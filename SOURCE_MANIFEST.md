# Source manifest

SHA-256 of every module under `src/`, as committed. `verify_snapshot.py`
exercises this code against recorded reference values; these digests let a
reader confirm the files have not changed since.

| file | sha256 |
|---|---|
| `src/lunaris/__init__.py` | `37cb9b4458186a7b42640345bec2987c54cf478c87b62fc44373ea330fb46c34` |
| `src/lunaris/_version.py` | `d1a7f7e129f2c0f0a3018ecf74fcf1584c0bff7a7695a724b84585ee42e3a80b` |
| `src/lunaris/common/__init__.py` | `3db3d816e7406c61705f7f1c0bd072208212b8a3db5203e80426fdb6fab9c4f9` |
| `src/lunaris/common/constants.py` | `cc3be2fed5f11218401ed728c7acaa048756510c57d11dafbfbfa7328e226ca4` |
| `src/lunaris/common/integrator_methods.py` | `3aa98b9c0847d1d709e1f0481c73fe6faba1e96d7680182f1cd662625b323ce0` |
| `src/lunaris/common/lunar_data.py` | `1807f0c3b3ebcab69ba48f2b381aafbc34ba322fd754ecff5bf2da89ac16be3d` |
| `src/lunaris/common/math_utils.py` | `55e69638bae3a734a4ea97f4ac30696ba7c8ed7161846f09455908c7406d6311` |
| `src/lunaris/common/paths.py` | `b43a340d5a3d95fcc15b9926dff0d03458caefb242a899a36d3b74fe1cb94f13` |
| `src/lunaris/common/provenance.py` | `59f4d175198870e92fc88fa2acb185b045c53402c702a505e0b4c2d5754db67d` |
| `src/lunaris/common/type_defs.py` | `95124205777313d706b29bc8ecbc90d451f1d499fc23cbdf6d0e55533a15c3c7` |
| `src/lunaris/core/__init__.py` | `65695ed22fa0dd8601b0d7df6b0a48572be3d184dd08af79aa310dc1a8d659f4` |
| `src/lunaris/core/propagation/__init__.py` | `65695ed22fa0dd8601b0d7df6b0a48572be3d184dd08af79aa310dc1a8d659f4` |
| `src/lunaris/core/propagation/integrators/__init__.py` | `63701c88a0e3b613868c4ee288da38850b4df420a90e6ed686e0fbb887fe4826` |
| `src/lunaris/core/propagation/integrators/symplectic.py` | `9196dc3485738c87546fc62ae2e8a03287111fed362e289d41d83658c45ced63` |
| `src/lunaris/loaders/__init__.py` | `bdbb5deb2cb0dee6402c3e205810a622bd6e936429a331bcdcdda15aea15fe5c` |
| `src/lunaris/loaders/io_gravity.py` | `716e034da0a102d2645ea1beffdf96c880c1c4597b1b468934c87a80de3d5677` |
| `src/lunaris/loaders/io_helpers.py` | `788a58e7c4fd9151e02824b34586e7b1466109e97fbaaae0ea90c07146e32e70` |
| `src/lunaris/loaders/spice_builder.py` | `438fc773245ad289247a5f915287294c5165e19489002629d37cf177526a8a15` |
| `src/lunaris/physics/__init__.py` | `68dc7b43990a3e3f7057d9af99c5d6226678acfe60dfd79f64ed2884c88de573` |
| `src/lunaris/physics/ephemeris.py` | `83030b6b74259dca574084c38015eb376a87714fa251ee35eb0552820716efa7` |
| `src/lunaris/physics/solar_effects.py` | `5004df62b22586aed8c382bee5cead8d3b81bd2542b5491e13bf3fde3088b0d7` |
| `src/lunaris/physics/solid_tides.py` | `3305dc368acfa1a47df573e0e53fb82cff5af0741de022406c466f5a0f219228` |
| `src/lunaris/physics/spherical_harmonics.py` | `6e2a0a85e18c9b400c003e6c2ad1e163c524c27f54bc0d404e0172709d16d740` |
| `src/lunaris/physics/third_body_effects.py` | `6665f42349c46e180c0a31ec139f9135bb4e72af0180fb5761ce9761cf7c2ed4` |
