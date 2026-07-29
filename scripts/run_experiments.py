"""Deney koşucusu: klasik SH çekirdeği makalesi için sayısal kanıt üretimi.

Tüm deneyler gerçek Lunaris deposu (D:\\Masaustu\\LUNAR_SIMULATION) ve gerçek
JGGRX 1800F katsayı dosyası üzerinde çalışır. Çıktılar bu makale klasörünün
metrics/ dizinine JSON olarak yazılır. Her deney tohumlanmıştır (seed) ve
kod yolu üretim çekirdeğinin kendisidir (GravityModel / sh_accel_* kernelleri).

Deneyler
--------
E1  Derece-bandı payları (band_share_analysis modülü ile, iki bant konfigürasyonu)
E2  Kesme (truncation) kriteri: tam spektrumdan ampirik N_min ve iki proxy modun
    (attenuation-only, Kaula p=1.7 / p=2.0) karşılaştırması
E3  Çekirdek zamanlaması: derece taraması + dual-pass vs 2x tek-geçiş
E4  Adaptif blend çekirdeği: C0 süreklilik taraması, curl (korunumsuzluk) ölçümü,
    ihmal edilen (U_hi - U_lo) grad(w) teriminin büyüklüğü
E5  Ayrık derece anahtarlama süreksizliği: kuantizasyon adımına göre ivme sıçraması

Kullanım (repo kökünden):
    python <bu_dosya> --out-dir <metrics_klasoru>
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(r"D:\Masaustu\LUNAR_SIMULATION")
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO))

from lunaris.common.lunar_data import resolve_lunar_gravity_path  # noqa: E402
from lunaris.common.math_utils import (  # noqa: E402
    LUNAR_DEGREE_POWER_EXPONENT,
    recommended_sh_degree,
)
from lunaris.physics.spherical_harmonics import (  # noqa: E402
    GravityModel,
    sh_accel_adaptive_blend_numba,
    sh_accel_fixed_numba,
    _compute_sh_acceleration_dual_numba,
    _apply_smoothstep,
)
from validation.gravity.band_share_analysis import (  # noqa: E402
    build_evidence_payload,
    compute_band_shares,
)

SEED = 20260719


def _dump(out_dir: Path, name: str, payload: dict) -> None:
    path = out_dir / name
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[yazildi] {path}")


# ---------------------------------------------------------------- E1: band payları
def exp1_band_shares(model_300: GravityModel, model_600: GravityModel,
                     gravity_file: Path, out_dir: Path) -> None:
    print("== E1: derece-bandi paylari ==")
    alts_main = [30.0, 50.0, 80.0, 100.0, 150.0, 200.0]
    results_main = [
        compute_band_shares(model_300, altitude_km=a, n_points=1000,
                            seed=SEED, band_edges=(60, 100))
        for a in alts_main
    ]
    payload_main = build_evidence_payload(
        gravity_file=gravity_file, model=model_300, band_edges=(60, 100),
        n_points=1000, seed=SEED, results=results_main,
    )
    _dump(out_dir, "e1_band_shares_60_100_nmax300.json", payload_main)

    # Kuyruk kontrolü: 100..300 bandı ve >300 kuyruğu, N_max=600 referansla
    alts_tail = [30.0, 80.0]
    results_tail = [
        compute_band_shares(model_600, altitude_km=a, n_points=400,
                            seed=SEED, band_edges=(100, 300))
        for a in alts_tail
    ]
    payload_tail = build_evidence_payload(
        gravity_file=gravity_file, model=model_600, band_edges=(100, 300),
        n_points=400, seed=SEED, results=results_tail,
    )
    _dump(out_dir, "e1_band_shares_100_300_nmax600.json", payload_tail)


# ------------------------------------------------------- E2: kesme kriteri karşılaştırması
def _degree_accel_rms(model: GravityModel, r_m: float) -> np.ndarray:
    """Gerçek katsayılardan derece başına ivme RMS'i (küre üzerinde, Kaula tipi).

    sigma_a(n) ~ (GM/r^2) (n+1) (R/r)^n sqrt(sum_m Cnm^2 + Snm^2)
    Tam-normalize katsayılar için küresel ortalama ivme katkısının standart
    büyüklük tahminidir; N_min karşılaştırmasında sadece oranlar kullanılır.
    """
    N = model.max_degree
    C = model.c_coeffs
    S = model.s_coeffs
    power = np.array([float(np.sum(C[n, : n + 1] ** 2 + S[n, : n + 1] ** 2))
                      for n in range(N + 1)])
    n_arr = np.arange(N + 1, dtype=np.float64)
    ratio = model.r_ref / r_m
    log_ratio_n = n_arr * math.log(ratio)
    amp = (model.mu / (r_m * r_m)) * (n_arr + 1.0) * np.exp(log_ratio_n) * np.sqrt(power)
    return amp


def _empirical_nmin(sigma: np.ndarray, tail_fraction: float) -> int:
    """Ampirik N_min: atılan kuyruk RMS payı esik altına ilk düştüğü N."""
    sq = sigma ** 2
    total = float(np.sum(sq[2:]))
    if total <= 0.0:
        return 0
    budget = (tail_fraction ** 2) * total
    tail = total
    for n in range(2, len(sq)):
        if tail <= budget:
            return n - 1
        tail -= sq[n]
    return len(sq) - 1


def exp2_truncation(model_1800: GravityModel, out_dir: Path) -> None:
    print("== E2: kesme kriteri ==")
    alts_km = [20.0, 30.0, 50.0, 80.0, 100.0, 150.0, 200.0, 300.0]
    rows = []
    for h in alts_km:
        r = model_1800.r_ref + h * 1000.0
        sigma = _degree_accel_rms(model_1800, r)
        n_emp = _empirical_nmin(sigma, 1e-2)
        n_att = recommended_sh_degree(h, model_1800.r_ref, attenuation_floor=1e-3)
        n_k17 = recommended_sh_degree(h, model_1800.r_ref,
                                      kaula_exponent=LUNAR_DEGREE_POWER_EXPONENT,
                                      kaula_tail_fraction=1e-2)
        n_k20 = recommended_sh_degree(h, model_1800.r_ref,
                                      kaula_exponent=2.0, kaula_tail_fraction=1e-2)
        rows.append({
            "altitude_km": h,
            "empirical_nmin_tail1e2": int(n_emp),
            "attenuation_only_floor1e3": int(n_att),
            "kaula_p1_7_tail1e2": int(n_k17),
            "kaula_p2_0_tail1e2": int(n_k20),
        })
        print(f"  h={h:6.0f} km  ampirik={n_emp:4d}  atten={n_att:4d}  "
              f"p1.7={n_k17:4d}  p2.0={n_k20:4d}")
    payload = {
        "seed": SEED,
        "model_max_degree": int(model_1800.max_degree),
        "tail_fraction": 1e-2,
        "attenuation_floor": 1e-3,
        "kaula_exponent_calibrated": LUNAR_DEGREE_POWER_EXPONENT,
        "note": "empirical_nmin gercek JGGRX katsayi spektrumundan; proxy modlar "
                "lunaris.common.math_utils.recommended_sh_degree",
        "rows": rows,
    }
    _dump(out_dir, "e2_truncation_criteria.json", payload)


# ---------------------------------------------------------------- E3: zamanlama
def _kernel_args(model: GravityModel):
    ws = model.make_workspace()
    return (model.r_ref, model.mu, model.c_coeffs, model.s_coeffs,
            model.diag_coeffs, model.subdiag_coeffs,
            model.a_coeffs, model.b_coeffs, model.scale_m_table,
            ws.P, ws.dP, ws.cos_m, ws.sin_m)


def exp3_timing(model_1800: GravityModel, out_dir: Path) -> None:
    print("== E3: zamanlama ==")
    args = _kernel_args(model_1800)
    r = model_1800.r_ref + 80e3
    lat, lon = math.radians(25.0), math.radians(40.0)
    x = r * math.cos(lat) * math.cos(lon)
    y = r * math.cos(lat) * math.sin(lon)
    z = r * math.sin(lat)

    degrees = [10, 20, 40, 60, 80, 120, 160, 200, 300, 400, 600, 800, 1200, 1800]
    # isinma (JIT)
    sh_accel_fixed_numba(x, y, z, 10, *args)

    rows = []
    for n in degrees:
        reps = max(5, min(400, int(4_000_000 / (n * n + 1))))
        times = []
        for _ in range(7):
            t0 = time.perf_counter_ns()
            for _ in range(reps):
                sh_accel_fixed_numba(x, y, z, n, *args)
            t1 = time.perf_counter_ns()
            times.append((t1 - t0) / reps)
        best_us = min(times) / 1000.0
        rows.append({"degree": n, "reps": reps, "best_us": best_us,
                     "median_us": statistics.median(times) / 1000.0})
        print(f"  N={n:5d}  {best_us:10.2f} us/eval")

    # dual-pass verimi: (60,120) tek geçiş vs iki ayrı çağrı
    _compute_sh_acceleration_dual_numba(x, y, z, 60, 120, *args)
    reps = 2000
    t0 = time.perf_counter_ns()
    for _ in range(reps):
        _compute_sh_acceleration_dual_numba(x, y, z, 60, 120, *args)
    dual_us = (time.perf_counter_ns() - t0) / reps / 1000.0
    t0 = time.perf_counter_ns()
    for _ in range(reps):
        sh_accel_fixed_numba(x, y, z, 60, *args)
        sh_accel_fixed_numba(x, y, z, 120, *args)
    two_us = (time.perf_counter_ns() - t0) / reps / 1000.0
    print(f"  dual(60,120)={dual_us:.2f} us   iki ayri cagri={two_us:.2f} us")

    payload = {
        "seed": SEED,
        "point": {"altitude_km": 80.0, "lat_deg": 25.0, "lon_deg": 40.0},
        "timer": "time.perf_counter_ns, best-of-7 blok, blok basina reps tekrar",
        "degree_sweep": rows,
        "dual_pass": {"pair": [60, 120], "dual_us": dual_us,
                      "two_single_calls_us": two_us,
                      "saving_pct": 100.0 * (1.0 - dual_us / two_us)},
    }
    _dump(out_dir, "e3_kernel_timing.json", payload)


# ------------------------------------------------- E4: blend çekirdeği analizi
def _blend_weight_chain(r: float, r_ref: float, n_min: int, n_max: int,
                        alt_far: float, alt_near: float, step: int):
    """Kernel'daki ağırlık zincirinin Python kopyası: (deg_lo, deg_hi, w)."""
    altitude = r - r_ref
    if altitude >= alt_far:
        return n_min, n_min, 0.0
    if altitude <= alt_near:
        return n_max, n_max, 0.0
    t = (alt_far - altitude) / (alt_far - alt_near)
    s = _apply_smoothstep(t)
    desired = n_min + s * (n_max - n_min)
    k = int(math.floor((desired - n_min) / step))
    deg_lo = min(max(n_min + k * step, n_min), n_max)
    deg_hi = min(deg_lo + step, n_max)
    if deg_hi == deg_lo:
        return deg_lo, deg_hi, 0.0
    w = min(max((desired - deg_lo) / (deg_hi - deg_lo), 0.0), 1.0)
    return deg_lo, deg_hi, w


def exp4_blend(model_300: GravityModel, out_dir: Path) -> None:
    print("== E4: adaptif blend analizi ==")
    m = model_300
    args = _kernel_args(m)
    n_far, n_near = 30, 120
    alt_far, alt_near = 200e3, 50e3
    step = 10
    lat, lon = math.radians(25.0), math.radians(40.0)
    u = np.array([math.cos(lat) * math.cos(lon),
                  math.cos(lat) * math.sin(lon),
                  math.sin(lat)])

    def blend_at(pos: np.ndarray) -> np.ndarray:
        ax, ay, az = sh_accel_adaptive_blend_numba(
            pos[0], pos[1], pos[2], n_far, n_near, alt_far, alt_near, step, *args)
        return np.array([ax, ay, az])

    def fixed_at(pos: np.ndarray, n: int) -> np.ndarray:
        ax, ay, az = sh_accel_fixed_numba(pos[0], pos[1], pos[2], n, *args)
        return np.array([ax, ay, az])

    # (a) radyal süreklilik taraması: 20..230 km, 4201 örnek (50 m adım)
    alts = np.linspace(20e3, 230e3, 4201)
    mags = np.empty_like(alts)
    for i, h in enumerate(alts):
        mags[i] = float(np.linalg.norm(blend_at((m.r_ref + h) * u)))
    jumps = np.abs(np.diff(mags))
    d_alt = float(alts[1] - alts[0])
    in_band = (alts[:-1] > alt_near) & (alts[1:] < alt_far)
    scan = {
        "d_alt_m": d_alt,
        "max_step_diff_in_band": float(np.max(jumps[in_band])),
        "max_step_diff_out_band": float(np.max(jumps[~in_band])),
        "note": "ardisik 50 m ornekler arasi |a| farki; C0 kontrolu",
    }
    print(f"  sureklilik: bant ici max adim farki {scan['max_step_diff_in_band']:.3e}, "
          f"bant disi {scan['max_step_diff_out_band']:.3e} m/s^2 (50 m adimda)")

    # (b) curl ölçümü: bant ortası noktada merkezi fark, h_fd = 0.5 m
    def curl_of(field, pos: np.ndarray, h_fd: float = 0.5) -> float:
        J = np.zeros((3, 3))
        for j in range(3):
            e = np.zeros(3); e[j] = h_fd
            J[:, j] = (field(pos + e) - field(pos - e)) / (2.0 * h_fd)
        c = np.array([J[2, 1] - J[1, 2], J[0, 2] - J[2, 0], J[1, 0] - J[0, 1]])
        return float(np.linalg.norm(c))

    pos_mid = (m.r_ref + 118e3) * u  # bant içi (rung ortasına yakın)
    curl_blend = curl_of(blend_at, pos_mid)
    curl_fixed = curl_of(lambda p: fixed_at(p, 120), pos_mid)
    a_pert_mid = float(np.linalg.norm(
        fixed_at(pos_mid, 120) + (m.mu / np.dot(pos_mid, pos_mid))
        * pos_mid / np.linalg.norm(pos_mid)))
    print(f"  curl: blend={curl_blend:.3e}  sabit-derece={curl_fixed:.3e}  [1/s^2]")

    # (c) ihmal edilen (U_hi - U_lo) * |dw/dr| terimi: bant boyunca tarama
    omitted_max = 0.0
    omitted_at_km = float("nan")
    for h in np.linspace(alt_near + 500.0, alt_far - 500.0, 300):
        r = m.r_ref + h
        dlo, dhi, w = _blend_weight_chain(r, m.r_ref, n_far, n_near,
                                          alt_far, alt_near, step)
        if dhi == dlo:
            continue
        dr = 25.0
        _, _, w_p = _blend_weight_chain(r + dr, m.r_ref, n_far, n_near,
                                        alt_far, alt_near, step)
        dlo_p, dhi_p, _ = _blend_weight_chain(r + dr, m.r_ref, n_far, n_near,
                                              alt_far, alt_near, step)
        if (dlo_p, dhi_p) != (dlo, dhi):
            continue  # basamak degisimi: w turevi tanimsiz, atla
        dw_dr = abs(w_p - w) / dr
        pos = r * u
        U_lo = m.potential_fixed(pos, degree=dlo)
        U_hi = m.potential_fixed(pos, degree=dhi)
        omitted = abs(U_hi - U_lo) * dw_dr
        if omitted > omitted_max:
            omitted_max = omitted
            omitted_at_km = h / 1000.0
    print(f"  ihmal edilen terim max |U_hi-U_lo||dw/dr| = {omitted_max:.3e} m/s^2 "
          f"(h={omitted_at_km:.1f} km); karsilastirma |a_pert|~{a_pert_mid:.3e}")

    payload = {
        "seed": SEED,
        "config": {"degree_far": n_far, "degree_near": n_near,
                   "alt_far_m": alt_far, "alt_near_m": alt_near,
                   "degree_step": step, "direction_lat_deg": 25.0,
                   "direction_lon_deg": 40.0},
        "continuity_scan": scan,
        "curl_1_s2": {"blend_mid_band": curl_blend,
                      "fixed_degree_reference": curl_fixed,
                      "fd_step_m": 0.5,
                      "point_altitude_km": 118.0},
        "omitted_term": {"max_m_s2": omitted_max, "at_altitude_km": omitted_at_km,
                         "pert_accel_scale_m_s2": a_pert_mid,
                         "dw_dr_fd_step_m": 25.0},
    }
    _dump(out_dir, "e4_blend_analysis.json", payload)


# ------------------------------------- E5: ayrık anahtarlama süreksizlik büyüklüğü
def exp5_switch_jump(model_300: GravityModel, out_dir: Path) -> None:
    print("== E5: ayrik anahtarlama sicramasi ==")
    m = model_300
    rng = np.random.default_rng(SEED)
    v = rng.normal(size=(400, 3))
    dirs = v / np.linalg.norm(v, axis=1, keepdims=True)
    ws = m.make_workspace()

    rows = []
    for h_km in [50.0, 100.0]:
        r = m.r_ref + h_km * 1000.0
        for base, q in [(120, 5), (120, 10), (120, 25), (120, 50)]:
            lo = base - q
            sq_jump = 0.0
            sq_pert = 0.0
            for u in dirs:
                pos = r * u
                a_hi = m.accel_fixed(pos, degree=base, workspace=ws)
                a_lo = m.accel_fixed(pos, degree=lo, workspace=ws)
                a_pm = -(m.mu / (r * r)) * u
                sq_jump += float(np.sum((a_hi - a_lo) ** 2))
                sq_pert += float(np.sum((a_hi - a_pm) ** 2))
            rms_jump = math.sqrt(sq_jump / len(dirs))
            rms_pert = math.sqrt(sq_pert / len(dirs))
            rows.append({"altitude_km": h_km, "degree_high": base, "step": q,
                         "rms_jump_m_s2": rms_jump,
                         "rms_pert_m_s2": rms_pert,
                         "jump_over_pert": rms_jump / rms_pert})
            print(f"  h={h_km:5.0f} km  {base}->{lo}: RMS sicrama {rms_jump:.3e} m/s^2 "
                  f"({100 * rms_jump / rms_pert:.3f}% pert)")
    payload = {"seed": SEED, "n_directions": 400, "rows": rows,
               "note": "RMS |a(N) - a(N-q)| yon ortalamasi; ayrik derece "
                       "anahtarlamasinin RHS'e enjekte ettigi sureksizlik olcegi"}
    _dump(out_dir, "e5_switch_jump.json", payload)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    gravity_file = resolve_lunar_gravity_path(None)
    print(f"gravity file: {gravity_file}")

    print("modeller yukleniyor (300 / 600 / 1800)...")
    model_300 = GravityModel.from_file(str(gravity_file), requested_degree=300)
    model_600 = GravityModel.from_file(str(gravity_file), requested_degree=600)
    model_1800 = GravityModel.from_file(str(gravity_file), requested_degree=1800)

    exp1_band_shares(model_300, model_600, Path(gravity_file), out_dir)
    exp2_truncation(model_1800, out_dir)
    exp3_timing(model_1800, out_dir)
    exp4_blend(model_300, out_dir)
    exp5_switch_jump(model_300, out_dir)
    print("tamamlandi.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
