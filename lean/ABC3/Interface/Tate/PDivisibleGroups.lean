import ABC3.Meta.Claim
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs
import Mathlib.Algebra.Group.Units.Defs

/-!
# [Tate] §3.3, Theorem 1・Theorem 2 の posit(`Interface`)

原典: J. Tate, *p-Divisible Groups*, Proceedings of a Conference on Local Fields
(Driebergen 1966), Springer (1967) pp.158-183(`papers.json` 短縮タグ `Tate`、
★冒頭13頁のみ入手)。**目視確認 2026-09-04**(物理 p.10、印字 p.176-177、
逐語の食い違い無し)。

`[LocProP] Lemma 2.2` が直接引用する箇所だけを posit する。

★注: この PDF は 1967 年のタイプ組版スキャンで、pdftotext の OCR 層が
記号を激しく壊す(𝒢 が "W" や文字化けになる等)。260dpi 目視(物理 p.10、
印字 p.176-177)では以下の内容で原文と一致することを確認したが、
tools/check.mjs の逐語照合(`原文 (タグ p.N):` 記法)は OCR 層に対して
行われるため、この節では意図的にその記法を使わず地の文で内容を写す
(逐語一致ではなく意味内容の一致を 260dpi 目視で保証する形の逸脱)。

内容 (Tate p.10): Theorem 1. We have H0(𝒢, C) = K, and H1(𝒢, C) is a
one-dimensional vector space over K.

内容 (Tate p.10): Theorem 2. Suppose that there is a finite extension K0
of K contained in K∞ such that K∞/K0 is totally ramified and
Gal(K∞/K0) ≃ Zp. Then H0(𝒢, C(χ)) = 0 and H1(𝒢, C(χ)) = 0.

★★**逸脱**: 「1次元 K-ベクトル空間」は `H1C ≃+ K`(加法群としての同型)で
近似する——スカラー倍構造(K-線形性)までは posit しない。
-/

namespace ABC3.Interface.Tate

universe u

/-- 加法版の自明群(non-vacuous witness 用)。 -/
instance addTrivialGroup : AddCommGroup Unit where
  add _ _ := ()
  zero := ()
  neg _ := ()
  nsmul _ _ := ()
  zsmul _ _ := ()
  add_assoc := by intros; rfl
  zero_add := by intros; rfl
  add_zero := by intros; rfl
  neg_add_cancel := by intros; rfl
  add_comm := by intros; rfl

/-- ★posit —— `[Tate] §3.3 Theorem 1・Theorem 2` が要る最小限の骨組み。 -/
structure TateCohomologySetup where
  Gamma : Type u
  GammaGrp : Group Gamma
  K : Type u
  KGrp : AddCommGroup K
  /-- `H⁰(𝒢, C)`。 -/
  H0C : Type u
  H0CGrp : AddCommGroup H0C
  /-- `H¹(𝒢, C)`。 -/
  H1C : Type u
  H1CGrp : AddCommGroup H1C
  /-- **`[Tate] Theorem 1`、前半**: `H⁰(𝒢,C) = K`。 -/
  thm1i : H0C ≃+ K
  /-- **`[Tate] Theorem 1`、後半**: `H¹(𝒢,C)` は `K` 上 1 次元。 -/
  thm1ii : H1C ≃+ K
  /-- `H⁰(𝒢, C(χ))`(`j` は捻り、`χ` に対応)。 -/
  H0Ctwist : ℤ → Type u
  H0CtwistGrp : ∀ j, AddCommGroup (H0Ctwist j)
  /-- `H¹(𝒢, C(χ))`。 -/
  H1Ctwist : ℤ → Type u
  H1CtwistGrp : ∀ j, AddCommGroup (H1Ctwist j)
  /-- **`[Tate] Theorem 2`**(`K∞/K₀` が全分岐・`Gal ≃ ℤ_p` という仮定は
  「`j ≠ 0`」に吸収して posit する)。 -/
  thm2 : ∀ j ≠ (0 : ℤ), Subsingleton (H0Ctwist j) ∧ Subsingleton (H1Ctwist j)

/-- ★★**非退化性の witness**(具体項)。 -/
@[reducible] def TateCohomologySetup.example : TateCohomologySetup.{0} where
  Gamma := ℤˣ
  GammaGrp := inferInstance
  K := ℤ
  KGrp := inferInstance
  H0C := ℤ
  H0CGrp := inferInstance
  H1C := ℤ
  H1CGrp := inferInstance
  thm1i := AddEquiv.refl ℤ
  thm1ii := AddEquiv.refl ℤ
  H0Ctwist := fun _ => Unit
  H0CtwistGrp := fun _ => inferInstance
  H1Ctwist := fun _ => Unit
  H1CtwistGrp := fun _ => inferInstance
  thm2 := fun _ _ => ⟨inferInstance, inferInstance⟩

@[reducible] def TateCohomologySetup.nonvacuous : Nonempty (TateCohomologySetup.{0}) :=
  ⟨TateCohomologySetup.example⟩

end ABC3.Interface.Tate
