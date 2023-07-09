(lang dune 3.9)

(using menhir 2.1)

(name min-caml)

(formatting
 (enabled_for dune))

(implicit_transitive_deps false)

(license BSD-3-Clause)

(maintainers "Eijiro Sumii" "Moe Masuko" "Kenichi Asai")

(authors "Eijiro Sumii" "Moe Masuko" "Kenichi Asai")

(source
 (github esumii/min-caml))

(generate_opam_files)

(package
 (name min-caml)
 (synopsis "An educational compiler for a minimal subset of OCaml")
 (description
  "An educational compiler for a minimal subset of OCaml, written in ~2000 lines of OCaml")
 (depends ocaml menhir))
