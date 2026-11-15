# Registration order

Both packages have new UUIDs, explicit dependency compatibility bounds, tests,
README files and their own copies of the MIT license. Find16S was an unregistered
prototype; its scratch data and original source are retained locally in ignored
files. This repository uses Julia's supported subdirectory-package layout.

1. Validate the intended commit in CI and inspect the real-genome comparison
   reports. Do not treat missing calls as harmless language-level differences.
2. Register ProfileHMMs first by commenting on that commit:
   `@JuliaRegistrator register subdir=packages/ProfileHMMs`
3. Wait for ProfileHMMs to merge into General. New packages normally have a
   three-day review window. This dependency must be registered before FindRNA
   can resolve through the registry.
4. Register the root package with `@JuliaRegistrator register` on the validated
   FindRNA commit, then check General's AutoMerge output.
5. TagBot creates `ProfileHMMs-v0.1.0` for the subpackage and `v0.1.0` for FindRNA
   after registration; do not manufacture a successful registration by tagging
   before the registry merge.

Registrator must have access to this repository. A commit comment alone is not
proof of registration: retain the resulting General pull-request URL and status.

References:
- https://github.com/JuliaRegistries/Registrator.jl#registering-a-package-in-a-subdirectory
- https://github.com/JuliaRegistries/General#automatic-merging-of-pull-requests
- https://github.com/JuliaRegistries/TagBot#subdirectory-packages
