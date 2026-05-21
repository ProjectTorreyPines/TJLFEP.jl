using Pkg
Pkg.activate("..")
Pkg.instantiate()
# using CUDA
# using TJLF
using FUSE
using IMAS
using GACODE
using IMASdd
using TurbulentTransport
using PackageCompiler

# create_sysimage(
#     [:FUSE, :IMAS, :TJLF, :GACODE, :IMASdd, :TurbulentTransport],
#     sysimage_path="TJLFEP_sysimage.so"
# )
create_sysimage(
    [:FUSE, :IMAS, :GACODE, :IMASdd, :TurbulentTransport],
    sysimage_path="noTJLF_TJLFEP_sysimage.so"
)
# create_sysimage(
#     [:FUSE, :IMAS],
#     sysimage_path="FUSE_IMAS_ONLY_TJLFEP_sysimage.so"
# )
