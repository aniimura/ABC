import ABC3.Meta.Claim
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs

/-!
# [BK] Lemma 3.8.1・Example 3.11 の posit(`Interface`)

原典: S. Bloch, K. Kato, *L-Functions and Tamagawa Numbers of Motives*,
The Grothendieck Festschrift I, Progress in Math. vol.86 (1990) pp.333-400
(`papers.json` 短縮タグ `BK`)。物理 p.23(印字 p.355、Lemma 3.8.1)・
物理 p.29(印字 p.361、Example 3.11)。**260 dpi 目視確認 2026-09-04**
(逐語の食い違い無し)。

★注: この PDF は Grothendieck Festschrift(1990)のスキャンで、pdftotext
の OCR 層が記号を激しく壊す(`B^+_dR` が "Bhn" 等に文字化け、`K̄` の上線
が落ちて "I<" や単なる "K" になる等)。260dpi 目視では以下の内容で
原文と一致することを確認したが、tools/check.mjs の逐語照合(Lean 側
`原文 (タグ p.N):` 記法・HTML `.verbatim`)は OCR 層に対して行われる
ため、この節は Tate と同じ運用にする: Lean コメントは地の文
(`内容 (BK p.N):`)で写し、HTML 側 `.verbatim` は OCR でも生き残る
短い見出し部分だけに絞る。

`[LocProP] Lemma 4.1` が直接引用する箇所だけを posit する:
- Lemma 3.8.1: `H^1(K, B^+_dR ⊗ V) → H^1(K, B_dR ⊗ V)` の単射性。
- Example 3.11, 式 (3.11.2): `H^1_e(K,T) = H^1_f(K,T) = H^1_g(K,T)`
  (de Rham 表現に対する3種の Selmer 条件の一致)。

de Rham 表現・結晶周期環そのものは posit しない(Falt1・Tate と同じ
流儀)。
-/

namespace ABC3.Interface.BK

universe u

structure DeRhamComparisonSetup where
  /-- `H^1(K, B^+_dR ⊗ V)`。 -/
  H1BdRPlus : Type u
  H1BdRPlusGrp : AddCommGroup H1BdRPlus
  /-- `H^1(K, B_dR ⊗ V)`。 -/
  H1BdR : Type u
  H1BdRGrp : AddCommGroup H1BdR
  /-- **[BK] Lemma 3.8.1** の射。 -/
  lem381 : H1BdRPlus →+ H1BdR
  /-- **[BK] Lemma 3.8.1**: 上の射は単射。 -/
  lem381_injective : Function.Injective lem381
  /-- `H^1_e(K,T)`(exponential Selmer 条件)。 -/
  H1e : Type u
  H1eGrp : AddCommGroup H1e
  /-- `H^1_f(K,T)`(finite Selmer 条件)。 -/
  H1f : Type u
  H1fGrp : AddCommGroup H1f
  /-- `H^1_g(K,T)`(de Rham Selmer 条件)。 -/
  H1g : Type u
  H1gGrp : AddCommGroup H1g
  /-- **[BK] Example 3.11**, 式 (3.11.2) の前半: `H^1_e ≅ H^1_f`。 -/
  ex311_ef : H1e ≃+ H1f
  /-- **[BK] Example 3.11**, 式 (3.11.2) の後半: `H^1_f ≅ H^1_g`。 -/
  ex311_fg : H1f ≃+ H1g

@[reducible] def DeRhamComparisonSetup.example : DeRhamComparisonSetup.{0} where
  H1BdRPlus := ℤ
  H1BdRPlusGrp := inferInstance
  H1BdR := ℤ
  H1BdRGrp := inferInstance
  lem381 := AddMonoidHom.id ℤ
  lem381_injective := fun _ _ h => h
  H1e := ℤ
  H1eGrp := inferInstance
  H1f := ℤ
  H1fGrp := inferInstance
  H1g := ℤ
  H1gGrp := inferInstance
  ex311_ef := AddEquiv.refl ℤ
  ex311_fg := AddEquiv.refl ℤ

@[reducible] def DeRhamComparisonSetup.nonvacuous : Nonempty (DeRhamComparisonSetup.{0}) :=
  ⟨DeRhamComparisonSetup.example⟩

end ABC3.Interface.BK
