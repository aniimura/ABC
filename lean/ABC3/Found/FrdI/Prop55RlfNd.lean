/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55Std

/-!
# [FrdI] Proposition 5.5, (iii) —— 係数拡大は non-dilating 性を保つ(還元)

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★`scModel_standardType`(`Prop55Std.lean`)は 3 つの入力を取る:
`hfsmff`(`𝒟` が FSMFF)・`hngl`(`𝒞^rlf` が not group-like)・
**`hnd`(係数拡大した `Φ` が non-dilating)**。
★前 2 つは閉じた(`hngl` は `scModel_not_isOfGroupLikeType`)。本ファイルは
**3 つめを単系の層の 3 条へ還元する**。

## ★還元の形

完全化の側(`MonoidOn.pfOn_isNonDilatingOn`)は

1. `M^pf` は sharp なので `isNonDilating_of_sharp` で **`M^pf` の上の含意**に降りる
2. `Pf.of` は準素元を保ち、`≼` を**保ち反射する**
3. `M` の non-dilating が効いて `α = id`

という 3 段だった。★係数拡大 `M ↦ S ⊗_ℕ M` でも**同じ 3 段**である。
本ファイルは 2 段目を**仮引数**(`hprim` / `hrefl`)にして、
1・3 段目を閉じる。

★★★**残るのは 2 段目だけ**である。`S = ℝ≥0` で `Φ` が perf-factorial なら
`Φ(A)` は素点ごとの座標を持つので、`≼` の反射は
**成分ごとの実数の比較を有理数の比較に戻す**だけで出る。

## ★★★測って分かったこと(2026-08-24) —— 座標がどの模型に乗っているか

我々は `Φ^rlf` の模型を **3 通り**持っている:

1. **テンソル模型** `ScT ℝ≥0 Φ = ℝ≥0 ⊗_ℕ Φ` —— `𝒞^rlf`(`scModel`)が使う
2. **因子模型** `Rlf ι ⊆ (Prime M → ℝ≥0)` —— `Definition 2.4` の定義そのもの
3. **錐模型** `rlfCone M ⊆ ℝ ⊗_ℤ M^gp` —— `Def24RlfCone` / `Def24RlfPerf`

★★**perf-factorial の座標(`factorMap` / `perfOrd`)は 2・3 の上にある**。
したがって上の 2 段目は、1 の上で座標を使うために何か橋が要る。

## ★★★★★訂正(2026-08-25)—— `rlf-agree` は要らなかった

本ファイルは当初ここに「**1 と 2(あるいは 3)の一致(節点 `rlf-agree`)を
経由しないと座標が使えない**」と書いていた。★**これは強すぎる見立てだった。**

要るのは**一致(同型)ではなく片道の準同型 1 本**である ——

```
scToFactor : ℝ≥0 ⊗_ℕ M →+ (Prime M → ℝ≥0)      (r ⊗ m ↦ r • factorMap ι (Pf.of m))
```

これは `TensorProduct.liftAddHom` で作れる。詳しくは `Prop55RlfRefl.lean`。
★**2 段目のうち `hrefl` はそこで閉じた**(`mprec_of_mprec_sc_primary`)。
残るのは `hprim`(準素元の保存)だけである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w u2 v2

section ScNd

variable {S : Type} [CommSemiring S]

/-- ★★★★★**係数拡大は non-dilating 性を保つ**(単系の層、還元した形)。

★1 段目(`isNonDilating_of_sharp`)と 3 段目(`M` の non-dilating)を閉じ、
2 段目(`toSc` が準素元を保ち `≼` を反射すること)を仮引数で受ける。 -/
theorem scMap_isNonDilating {M : Type w} [AddCommMonoid M]
    (hs : IsSharp M) (hsS : IsSharp (ScT S M))
    (hprim : ∀ b : M, IsPrimaryElt b → ∃ a : ScT S M, IsPrimaryElt a
      ∧ MPrec (toSc (S := S) b) a ∧ MPrec a (toSc (S := S) b))
    (hrefl : ∀ b : M, IsPrimaryElt b → ∀ x : M,
      MPrec (toSc (S := S) x) (toSc b) → MPrec x b)
    (f : M →+ M) (h : IsNonDilating f) :
    IsNonDilating (scMap (S := S) f) := by
  refine isNonDilating_of_sharp hsS _ (fun key => ?_)
  -- ★`M` の側の含意を作る
  have hM : ∀ a : MChar M, IsPrimaryElt a → MPrec (charMap f a) a := by
    intro a ha
    set b : M := (mCharEquivOfSharp hs).symm a with hb
    have hab : toChar b = a := (mCharEquivOfSharp hs).apply_symm_apply a
    have hbp : IsPrimaryElt b := by
      have := isPrimaryElt_map (mCharEquivOfSharp hs).symm ha
      rwa [← hb] at this
    -- ★★`toSc b` そのものが準素元でなくてもよい —— **`≼`-同値な準素元があれば足りる**
    obtain ⟨z, hzp, h1, h2⟩ := hprim b hbp
    have hkey : MPrec (scMap (S := S) f (toSc b)) (toSc b) :=
      mprec_trans (mprec_trans (mprec_map (scMap (S := S) f) h1) (key z hzp)) h2
    rw [scMap_toSc] at hkey
    have hMb : MPrec (f b) b := hrefl b hbp _ hkey
    have := mprec_map (toChar : M →+ MChar M) hMb
    rw [← charMap_toChar f b, hab] at this
    exact this
  -- ★`M` の non-dilating から `charMap f = id`、sharp から `f = id`
  have hchar : charMap f = AddMonoidHom.id (MChar M) := h hM
  have hf : f = AddMonoidHom.id M := by
    ext x
    have h1 : toChar (f x) = toChar x := by
      have := congrArg (fun t : MChar M →+ MChar M => t (toChar x)) hchar
      rwa [charMap_toChar] at this
    exact toChar_injective_of_isSharp hs h1
  rw [hf, scMap_id]

/-- ★**強い `hprim`(「`toSc b` が準素元」)は弱い形を含む**。

★★弱い形で十分なのは、`≼`-同値な準素元 `z` があれば
`scMap f (toSc b) ≼ scMap f z ≼ z ≼ toSc b` と繋がるからである。 -/
theorem hprimWeak_of_isPrimaryElt {M : Type w} [AddCommMonoid M]
    (hp : ∀ b : M, IsPrimaryElt b → IsPrimaryElt (toSc (S := S) b)) :
    ∀ b : M, IsPrimaryElt b → ∃ z : ScT S M, IsPrimaryElt z
      ∧ MPrec (toSc (S := S) b) z ∧ MPrec z (toSc (S := S) b) :=
  fun b hb => ⟨toSc b, hp b hb, mprec_refl _, mprec_refl _⟩

end ScNd

/-! ## ★2. `MonoidOn` の層 -/

section ScNdOn

variable {D : Type u} [Category.{v} D] {S : Type} [CommSemiring S]

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** ——
**係数拡大は non-dilating 性を保つ**(`MonoidOn` の層、還元した形)。

★これが `scModel_standardType` の残る入力 `hnd` である。 -/
theorem MonoidOn.scOn_isNonDilatingOn (Φ : MonoidOn.{v, u, w} D)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hs : ∀ A : D, IsSharp (Φ.val A)) (hsS : ∀ A : D, IsSharp (ScT S (Φ.val A)))
    (hprim : ∀ (A : D) (b : Φ.val A), IsPrimaryElt b →
      ∃ z : ScT S (Φ.val A), IsPrimaryElt z
        ∧ MPrec (toSc (S := S) b) z ∧ MPrec z (toSc (S := S) b))
    (hrefl : ∀ (A : D) (b : Φ.val A), IsPrimaryElt b → ∀ x : Φ.val A,
      MPrec (toSc (S := S) x) (toSc b) → MPrec x b)
    (h : Φ.IsNonDilatingOn) : (phiScOn S Φ hcharInj).IsNonDilatingOn := by
  intro A e
  show IsNonDilating (scMap (S := S) (Φ.map e))
  exact scMap_isNonDilating (hs A) (hsS A) (hprim A) (hrefl A) _ (h A e)

end ScNdOn

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Proposition 5.5, (iii)` の「係数拡大は non-dilating 性を保つ」
(単系の層の 2 条 —— 準素元の保存と `≼` の反射 —— へ還元した形)。 -/
def MonoidOn.scOn_isNonDilatingOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 係数拡大は non-dilating 性を保つ(単系の 2 条へ還元)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
