# How to release varr-numerics

* increment the version number in `CMakeLists.txt`
* add changelog entries for the new version in
  * `debian/changelog`
  * `varr-numerics.spec.in`
* create a tarball

  ```shell
  cmake .
  cmake --build . --target package_source
  ```

  this will create `varr-numerics-{version}.tar.xz`
