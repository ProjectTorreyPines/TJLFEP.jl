#!/bin/bash -l
# Stage the full TJLFEP build tree onto CFS so sysimages can be baked with (and run from)
# m3739-readable paths only. Run on a login node as the maintainer:
#
#   bash build/sysimage/publish_build_tree.sh
#
# Stages the dev repos the FUSE manifest path-deps on (working trees, INCLUDING
# uncommitted mods and .git/ for provenance) into $CFS_DIR/src/<name>, rewrites the
# absolute path deps in the staged FUSE Manifest to relative sibling paths, and records
# per-repo git state in src/BUILD_STATE. Idempotent: re-run to refresh (rsync --delete).
# IMASdd/IMASggd resolve from FuseRegistry (freed 2026-08-25) and are not staged.
#
# Perms: the original sharing bug was `rsync -a`/`cp` preserving the maintainer's primary
# group. Here rsync runs WITHOUT -g/-o so setgid parent dirs assign group m3739, with an
# explicit chgrp/setgid pass as belt-and-braces. Staged tree is group read-only (g+rX);
# the maintainer keeps owner write for rebuilds.

set -euo pipefail
umask 007

CFS_DIR="${TJLFEP_CFS_DIR:-/global/cfs/cdirs/m3739/TJLFEP}"
SRC="${CFS_DIR}/src"
SCRATCH_DEV="${SCRATCH_DEV:-/pscratch/sd/t/tneiser/.julia/dev}"
HOME_DEV="${HOME_DEV:-${HOME}/.julia/dev}"
GROUP=m3739

SCRATCH_REPOS=(TJLFEP TJLF FUSE ALPHA TurbulentTransport)
HOME_REPOS=(CHEASE TroyonBetaNN)

mkdir -p "${SRC}"

stage() {
    local from="$1" name="$2"
    echo "--- staging ${name} from ${from}/${name}"
    rsync -rlpt --delete \
        --chmod=Dg+s,g-w,g+rX,o-rwx \
        --exclude='*.so' --exclude='*.so.bak' --exclude='*.out' --exclude='*.err' \
        --exclude='smoke_out_*/' --exclude='tjlfep_smoke_out_*/' --exclude='attic/' \
        --exclude='build/ad/' --exclude='ucp_fileInput_*/' \
        --exclude='.git/objects/pack/tmp_*' \
        "${from}/${name}/" "${SRC}/${name}/"
}

for r in "${SCRATCH_REPOS[@]}"; do stage "${SCRATCH_DEV}" "$r"; done
for r in "${HOME_REPOS[@]}";    do stage "${HOME_DEV}"    "$r"; done

# Rewrite the FUSE manifest's absolute dev paths (maintainer scratch/home) to relative
# sibling paths, which resolve against the staged project dir for every user.
FMAN="${SRC}/FUSE/Manifest.toml"
for r in CHEASE TroyonBetaNN; do
    sed -i "s|^path = \".*/\.julia/dev/${r}\"$|path = \"../${r}\"|" "${FMAN}"
done
if grep -n 'path = "/' "${FMAN}"; then
    echo "ERROR: absolute path deps remain in ${FMAN}" >&2
    exit 1
fi
echo "--- FUSE manifest path deps:"
grep '^path = ' "${FMAN}"

# Record staged provenance.
{
    echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    for r in "${SCRATCH_REPOS[@]}"; do
        echo "${r}: $(git -C "${SCRATCH_DEV}/${r}" describe --always --dirty --tags 2>/dev/null || echo unknown)"
    done
    for r in "${HOME_REPOS[@]}"; do
        echo "${r}: $(git -C "${HOME_DEV}/${r}" describe --always --dirty --tags 2>/dev/null || echo unknown)"
    done
} > "${SRC}/BUILD_STATE"
cat "${SRC}/BUILD_STATE"

# Belt-and-braces group repair (rsync without -g relies on setgid inheritance).
chgrp -R "${GROUP}" "${SRC}"
find "${SRC}" -type d ! -perm -g+s -exec chmod g+s {} +

BAD=$(find "${CFS_DIR}" ! -group "${GROUP}" | head -5)
if [[ -n "${BAD}" ]]; then
    echo "ERROR: entries not group ${GROUP}:" >&2
    echo "${BAD}" >&2
    exit 1
fi
echo "=== build tree staged at ${SRC} (all group ${GROUP}) ==="
