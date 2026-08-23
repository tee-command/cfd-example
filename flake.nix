{
  description = "C++ development environment (cmake, clang, eigen, etc.)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            # Development Tools
            cmake
            gnumake
            clang-tools
            clang
            eigen
            git
            ninja
            llvmPackages.openmp
            pkg-config
          ];
        };
      });
}
