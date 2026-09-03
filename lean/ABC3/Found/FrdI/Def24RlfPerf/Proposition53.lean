/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24RlfCone
import ABC3.Found.FrdI.Def24
import ABC3.Found.FrdI.Thm42
import ABC3.Found.FrdI.Def24RlfPerf.Definition24

/-!
# Def24RlfPerf —— `[FrdI] Proposition 5.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct NNReal
universe v u w
variable {M : Type w} [AddCommMonoid M]

/-! ## ★5. `Φ^rlf : MonoidOn 𝒟` —— 仮定なし -/

variable {D : Type u} [Category.{v} D]

/-- ★★★★★★★**perf-factorial な `Φ` の実化は仮定なしで `MonoidOn 𝒟` になる**。

★`Proposition 5.3` は `Φ` が perf-factorial であることを**仮定している**ので、
これで `Φ^rlf` は完全に立つ。 -/
noncomputable def phiRlfConeOnOfPerfFactorial (Φ : MonoidOn.{v, u, w} D)
    (hdiv : Φ.IsDivisorialOn) (hpf : Φ.IsPerfFactorialOn) : MonoidOn.{v, u, w} D :=
  phiRlfConeOnOfDivisorial Φ hdiv (fun A => isSharp_rlfCone_of_perfFactorial (hpf A))

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Definition 2.4, (i)` の実化が perf-factorial な `Φ` について
`MonoidOn 𝒟` になること(`Proposition 5.3` が要求するもの)。 -/
def phiRlfConeOnOfPerfFactorial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — Φ^rlf(perf-factorial なら仮定なしで MonoidOn 𝒟)",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
