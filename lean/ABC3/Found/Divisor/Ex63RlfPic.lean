/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex63RlfRatStd
import ABC3.Found.Divisor.ArithPicR

/-!
# [FrdI] Theorem 6.4, (i) —— `δ_A : Pic_Φ(A) ≅ ℝ` を `Φ^rlf` の水準で述べる

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> which is an immediate consequence of the well-known Dirichlet unit theorem

## ★★何をする節点か

`ArithPicR.lean` は **Dirichlet 単数定理 ＋ 類数有限**から

    principalSpan L = LinearMap.ker arithDegreeLin
    arithPicIso : ((ArithPlace L →₀ ℝ) ⧸ principalSpan L) ≃ₗ[ℝ] ℝ

を出した（節点 `thm64-i-dirichlet`）。★本ファイルはこれを
**`𝒞^rlf` の `Φ^gp` の水準**に載せ替える。

## ★★★`Φ^rlf(L)^gp` は「全実係数の因子」そのもの

`Ex63Rlf.lean` の `arithDatumRlf`（`grp = ⊤`）では

    Φ^rlf(L) = effR ⊤ = ℝ≥0[V(L)]  ⟹  Φ^rlf(L)^gp ≃+ ℝ[V(L)]

（`effRGpEquiv` ＋ `AddSubgroup.topEquiv`）。★★したがって `arithPicIso` の定義域と
**ちょうど噛み合う**。

## ★★★★`Φ^birat` は**実現化した**主因子である

`Pic_Φ(A) = Φ^gp(A)/Φ^birat(A)` で、`𝒞^rlf` では `Φ^birat` も実現化されるので
**`principalSpan L`（主因子の ℝ-span）**である。
★`Ex63Rlf.lean` の `rlfModelData` は有理関数の単系として `B(L) = L^×`（ℤ 係数）を
そのまま使っているので、そこから来る birational 部分は ℤ-span である。
★★本ファイルの `phiBiratRlf` は**実現化した側**（ℝ-span）を明示的に置く ——
これが原文の `𝒞^rlf` における `Φ^birat` である。

★★★★★述べ方は**第一同型定理の形**にした（商の型の食い違いを避けるため）:
`rlfDeg` が全射で、その核が `phiBiratRlf` である。`rlfDeltaA` がそこから出る同型。
-/

namespace ABC3.Found.Divisor

open ABC3.Found.FrdI CategoryTheory _root_.NumberField
open ABC3.Meta

section

variable (F Kbar : Type) [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
  [IsGalois F Kbar]

/-! ## ★1. `Φ^rlf(L)^gp ≃+ ℝ[V(L)]` -/

/-- ★★**`Φ^rlf(L)^gp ≃+ ℝ[V(L)]`**（全実係数の因子）。

★`Φ^rlf(L) = effR ⊤` なので、`effRGpEquiv` で `⊤` に落ち、`topEquiv` で全体になる。 -/
noncomputable def rlfGpEquiv (Y : (FinSub F Kbar)ᵒᵖ) :
    Gp (((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).val Y)
      ≃+ (ArithPlace Y.unop.toIF →₀ ℝ) :=
  (effRGpEquiv (⊤ : AddSubgroup (ArithPlace Y.unop.toIF →₀ ℝ))
    isGenSubgroupR_top).trans (AddSubgroup.topEquiv)

/-! ## ★2. `Φ^birat`（実現化した主因子） -/

/-- ★★**`Φ^rlf` の birational 部分** —— 実現化した主因子の張る部分群
（`principalSpan L` を `rlfGpEquiv` で引き戻したもの）。 -/
noncomputable def phiBiratRlf (Y : (FinSub F Kbar)ᵒᵖ) :
    AddSubgroup (Gp (((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).val Y)) :=
  AddSubgroup.comap (rlfGpEquiv F Kbar Y).toAddMonoidHom
    (principalSpan Y.unop.toIF).toAddSubgroup

/-! ## ★3. `δ_A` -/

/-- ★★**`δ_A : Φ^rlf(A)^gp → ℝ`** —— 算術因子の次数。 -/
noncomputable def rlfDeg (Y : (FinSub F Kbar)ᵒᵖ) :
    Gp (((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).val Y) →+ ℝ :=
  (arithDegreeLin (L := Y.unop.toIF)).toAddMonoidHom.comp (rlfGpEquiv F Kbar Y).toAddMonoidHom

/-- ★**`δ_A` は全射**（無限素点に台をもつ因子で任意の値を取れる）。 -/
theorem rlfDeg_surjective (Y : (FinSub F Kbar)ᵒᵖ) :
    Function.Surjective (rlfDeg F Kbar Y) := by
  intro c
  obtain ⟨x, hx⟩ := arithDegreeLin_surjective (L := Y.unop.toIF) c
  exact ⟨(rlfGpEquiv F Kbar Y).symm x, by
    show arithDegreeLin ((rlfGpEquiv F Kbar Y) ((rlfGpEquiv F Kbar Y).symm x)) = c
    rw [AddEquiv.apply_symm_apply]
    exact hx⟩

/-- ★★★★★★**`δ_A` の核はちょうど `Φ^birat`** ——
これが **Dirichlet 単数定理 ＋ 類数有限**の中身（`principalSpan_eq_ker`）である。 -/
theorem rlfDeg_ker (Y : (FinSub F Kbar)ᵒᵖ) :
    (rlfDeg F Kbar Y).ker = phiBiratRlf F Kbar Y := by
  ext z
  constructor
  · intro hz
    show (rlfGpEquiv F Kbar Y) z ∈ (principalSpan Y.unop.toIF).toAddSubgroup
    rw [principalSpan_eq_ker]
    exact hz
  · intro hz
    have h1 : (rlfGpEquiv F Kbar Y) z ∈ (principalSpan Y.unop.toIF).toAddSubgroup := hz
    rw [principalSpan_eq_ker] at h1
    exact h1

/-- ★★★★★★★**[FrdI] Theorem 6.4, (i)** —— `δ_A : Pic_Φ(A) ≃ ℝ`。

原文 (FrdI p.114):
> which is an immediate consequence of the well-known Dirichlet unit theorem

★★中身は `rlfDeg` の全射性（無限素点）と、核が `Φ^birat` であること
（**Dirichlet 単数定理 ＋ 類数有限**、`ArithPicR.lean`）だけである。 -/
noncomputable def rlfDeltaA (Y : (FinSub F Kbar)ᵒᵖ) :
    (Gp (((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType).val Y) ⧸ phiBiratRlf F Kbar Y)
      ≃+ ℝ :=
  (QuotientAddGroup.quotientAddEquivOfEq (rlfDeg_ker F Kbar Y).symm).trans
    (QuotientAddGroup.quotientKerEquivOfSurjective _ (rlfDeg_surjective F Kbar Y))

end

/-! ### ★出典の紐付け -/

def rlfDeltaA.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — δ_A : Pic_Φ(A) ≃ ℝ(Φ^rlf の水準)",
    sectionId := "frdi-thm-6-4" }

def rlfDeltaA.needs : List ProofObligation :=
  [ .citation "[ABC3]" "principalSpan_eq_ker(Dirichlet 単数定理 ＋ 類数有限)"
      (.inProject "ABC3" "ABC3.Found.Divisor.principalSpan_eq_ker") 114,
    .citation "[ABC3]" "arithDegreeLin_surjective(無限素点で任意の次数が取れる)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDegreeLin_surjective") 114,
    .implicitStep
      ("★Φ^birat は実現化した主因子(ℝ-span)である。rlfModelData の有理関数の単系は " ++
       "B(L) = L^×(ℤ 係数)なので、そこから来る birational 部分は ℤ-span であり、" ++
       "本ファイルの phiBiratRlf はその実現化を明示的に置いたものである") 114 ]

def rlfGpEquiv.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Φ^rlf(L)^gp ≃+ ℝ[V(L)]",
    sectionId := "frdi-thm-6-4" }

def phiBiratRlf.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Φ^rlf の birational 部分(実現化した主因子)",
    sectionId := "frdi-thm-6-4" }

def rlfDeg.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — δ_A : Φ^rlf(A)^gp → ℝ(算術因子の次数)",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.Divisor
