/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55Pf
import ABC3.Found.FrdI.Prop32Equiv
import ABC3.Found.FrdI.Prop53UntrPf
import Mathlib.CategoryTheory.Conj

/-!
# [FrdI] Proposition 5.5, (i) の**零因子**の側 —— `𝒞^pf` の `Div` は分母を持つ

原文 (FrdI p.104):
> Proof. Assertion (i) follows immediately for Frobenius-trivial A by considering

★★★**なぜこれが要るか**: `Proposition 5.3` は `(𝒞^un-tr)^pf` の有理関数の単系を
`ℚ·Φ^birat = Φ^birat ⊗_ℤ ℚ = (Φ^birat)^pf` と書き、`Proposition 5.5, (iv)` は
それを「definitions から immediate」と畳む。★我々の側では `𝒞^pf` の有理関数の単系は
`(Φ^pf)^birat` として作られているので、**`ℚ` がどこから出るのか**を明示する必要がある。

★★★★**答: `𝒞^pf` の零因子の分母から出る**。
「`α` の `k` 乗根」にあたる `𝒞^pf` の自己射の零因子は `Div α / k` である
(`rootDiv_otriPfMk`)。★これが本ファイルの中身であり、
`Div(𝒪^▷(A^pf)) = ℚ≥0 · Div(𝒪^▷(A))`(`rootDiv_otri_image`)がその帰結である。

## ★筋(3 本)

1. `Div_conjRt_pull` —— 共役 `conjRt` は零因子を動かさない
   (`rtExt A 1` は**同型**なので `Div = 0`、次数 `1`)。
2. `rootDiv_zetaMk` / `rootDiv_otriPfMk` —— 代表元の零因子の明示式。
   ★`rootDiv_mk_sameRoot`(同根の `Div` の明示式)に、
   添字 `ζ` の**底恒等性**と**次数**を入れるだけである。
3. `rootDiv_otri_image` —— 上の 2 本を、在庫の
   `otri_pfRoot_exists_rep`(`𝒪^▷(A^pf)` の元はすべて `k` 乗根)と
   `otriPfMk_mem`(`k` 乗根は `𝒪^▷` に入る)で両側から挟む。

## ★★逸脱・限界(記録)

* 本ファイルは **`𝒪^▷` の層**の主張である。原文が要求するのは `Φ^birat`
  (すなわち `𝒪^×(A^birat)` の `Div^gp` の像)なので、
  ここから `𝒞^birat` へ移す段(`P := biratPre P G` で読み直す)が残る。
* `𝒪^▷` から `𝒪^×` への段、および `Div` から `Div^gp` への段も残る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★1. 共役は零因子を動かさない -/

/-- ★★**共役 `conjRt` は零因子を動かさない**(`rtExt A 1` へ引き戻せば元に戻る)。

★`rtExt A 1` は**同型**(`isIso_rtExt_one`)なので `Div = 0`・次数 `1` であり、
`Div_comp` を 2 回使うと `Φ.map (Base (inv (rtExt A 1)))` だけが残る。 -/
theorem Div_conjRt_pull {A : C} (f : A ⟶ A) :
    (Φ.map (P.Base (rtExt P F A 1))) (P.Div (conjRt (F := F) f)) = P.Div f := by
  haveI := isIso_rtExt_one P F A
  have hA0 : P.Div (rtExt P F A 1) = 0 := (rtExt_frobType P F A 1).1.2
  have hinv0 : P.Div (inv (rtExt P F A 1)) = 0 := Div_inv_eq_zero _ hA0
  have hx : P.Div (conjRt (F := F) f)
      = Φ.map (P.Base (inv (rtExt P F A 1))) (P.Div f) := by
    show P.Div (inv (rtExt P F A 1) ≫ f ≫ rtExt P F A 1) = _
    rw [P.Div_comp, P.Div_comp, hA0, hinv0]
    simp [rtExt_degFr]
  rw [hx, ← Φ.map_comp (P.Base (inv (rtExt P F A 1))) (P.Base (rtExt P F A 1)) (P.Div f),
    ← P.Base_comp, IsIso.hom_inv_id, P.Base_id, Φ.map_id]

/-! ## ★2. 代表元の零因子の明示式 —— **分母が出る** -/

/-- ★★★★★**`ζ` 添字の代表元の零因子** —— 分母は**添字の Frobenius 次数**。

★これが「`𝒞^pf` の零因子は `Φ` の元を `ℕ≥1` で割ったもの」という事実の本体である。 -/
theorem rootDiv_zetaMk {A : C} {z : rtObj P F A 1 ⟶ rtObj P F A 1}
    (hz : IsFrobeniusType P z) (hbz : IsBaseIdentity P z)
    (φ : rtObj P F A 1 ⟶ rtObj P F A 1) :
    rootDiv (zetaMk (F := F) hz φ)
      = Pf.mk (Φ.map (P.Base (rtExt P F A 1)) (P.Div φ)) (P.degFr z) := by
  have h := rootDiv_mk_sameRoot (P := P) (F := F) (r := 1) (idxZeta' (F := F) hz) φ
  have h2 : (Φ.map (P.Base z)) (P.Div φ) = P.Div φ := by
    rw [show P.Base z = P.Base (𝟙 _) from hbz, P.Base_id, Φ.map_id]
  refine (show rootDiv (zetaMk (F := F) hz φ)
      = (pfRootPre P F).Div (show HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨A, 1⟩ from
        HomPf.mk (idxZeta' (F := F) hz) φ) from rfl).trans (h.trans ?_)
  refine congrArg₂ Pf.mk ?_ ?_
  · exact congrArg (fun t : Φ.val (P.toElem.obj (rtObj P F A 1)).base =>
      (Φ.map (P.Base (rtExt P F A 1))) t) h2
  · show 1 * 1 * repRoot (idxZeta' (F := F) hz) = P.degFr z
    rw [one_mul, one_mul]
    rfl

/-- ★★★★★**「`α` の `k` 乗根」の零因子は `Div α / k`**。

★★**ここが `ℚ` の出どころ**である —— `𝒞` の側には無い「`k` で割る」が、
`𝒞^pf` では添字の Frobenius 次数として現れる。 -/
theorem rootDiv_otriPfMk (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ζ : A ⟶ A} (hζ : IsFrobeniusType P ζ) (hbz : IsBaseIdentity P ζ) (α : End A) :
    rootDiv (otriPfMk (F := F) hiso hζ α) = Pf.mk (P.Div (α : A ⟶ A)) (P.degFr ζ) := by
  have key : (Φ.map (P.Base (rtExt P F A 1)))
      ((Φ.map (P.Base (conjRt (F := F) ζ))) (P.Div (conjRt (F := F) ((α : A ⟶ A)))))
      = P.Div (α : A ⟶ A) := by
    have h1 : P.Base (conjRt (F := F) ζ) = 𝟙 (P.toElem.obj (rtObj P F A 1)).base := by
      have h0 := conjRt_baseIdentity (F := F) hbz
      rw [show P.Base (conjRt (F := F) ((ζ : A ⟶ A))) = P.Base (𝟙 _) from h0, P.Base_id]
    rw [h1, Φ.map_id]
    exact Div_conjRt_pull _
  have h := rootDiv_mk_sameRoot (P := P) (F := F) (r := 1)
    (idxZeta hiso hζ) (conjRt (F := F) ((α : A ⟶ A)))
  refine (show rootDiv (otriPfMk (F := F) hiso hζ α)
      = (pfRootPre P F).Div (show HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨A, 1⟩ from
        HomPf.mk (idxZeta hiso hζ) (conjRt (F := F) ((α : A ⟶ A)))) from rfl).trans
    (h.trans ?_)
  refine congrArg₂ Pf.mk key ?_
  show 1 * 1 * repRoot (idxZeta (F := F) hiso hζ) = P.degFr ζ
  rw [one_mul, one_mul]
  exact conjRt_degFr ζ

/-! ## ★3. `Div(𝒪^▷(A^pf)) = ℚ≥0 · Div(𝒪^▷(A))` -/

/-- ★★★★★★**`𝒞^pf` の `𝒪^▷` の零因子の像は、`𝒞` の側の像の `ℚ≥0`-スパン**。

★★★これが `Proposition 5.3` の「`ℚ·Φ^birat`」の `ℚ` の出どころである
(まだ `𝒪^▷` の層。`Φ^birat` へ移すには `𝒞^birat` で読み直す段が要る)。

* `⊆` —— 在庫 `otri_pfRoot_exists_rep`(`𝒪^▷(A^pf)` の元はすべて `k` 乗根)＋
  `rootDiv_zetaMk`＋`Div_conjRt_pull`。
* `⊇` —— `otriPfMk_mem`(`k` 乗根は `𝒪^▷` に入る)＋ `rootDiv_otriPfMk`。 -/
theorem rootDiv_otri_image (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A)
    (ζ : ℕ+ →* End A) (hdeg : ∀ n : ℕ+, P.degFr ((ζ n : End A) : A ⟶ A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P ((ζ n : End A) : A ⟶ A)) :
    {d : Pf (Φ.val (P.toElem.obj A).base) |
        ∃ f : End (⟨A, 1⟩ : PfRootObj P F),
          f ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F) ∧ rootDiv f = d}
      = {d : Pf (Φ.val (P.toElem.obj A).base) |
        ∃ (α : End A) (k : ℕ+), α ∈ OTri P A ∧ d = Pf.mk (P.Div (α : A ⟶ A)) k} := by
  ext d
  constructor
  · rintro ⟨f, hf, rfl⟩
    obtain ⟨z, hz, hbz, φ, hφ, hmk⟩ := otri_pfRoot_exists_rep (F := F) hiso hA hf
    refine ⟨conjRtInv (F := F) φ, P.degFr z, conjRtInv_otri hφ, ?_⟩
    rw [← hmk, rootDiv_zetaMk hz hbz φ]
    refine congrArg₂ Pf.mk ?_ rfl
    rw [← conjRt_conjRtInv (F := F) φ, Div_conjRt_pull, conjRtInv_conjRt]
  · rintro ⟨α, k, hα, rfl⟩
    refine ⟨otriPfMk (F := F) hiso (hprop k).2 α, otriPfMk_mem hiso (hprop k).2 α hα, ?_⟩
    rw [rootDiv_otriPfMk hiso (hprop k).2 (hprop k).1 α, hdeg k]

/-! ## ★4. 分母を払う —— `𝒪^▷(A^pf)` の元は正の冪で `𝒞` へ戻る

★★★これが `(Φ^pf)^birat ⊆ ℚ·Φ^birat`(逆向きの包含)の**単系の側の骨**である。
`𝒪^▷(A^pf) ≅ 𝒪^▷(A)^pf`(`otriPfEquiv`、`Proposition 5.5, (i)`)の下で
`k · (α/k) = α`(`Pf.nsmul_mk_self`)を読み替えただけである。 -/

/-- ★★★★★**`𝒪^▷(A^pf)` の元 `f` には正の `k` があって `f^k` は `𝒞` の元の像**。

★`otriPfMap` の全射性で `f = α^{1/k}` と書き、
同型の加法性(`map_nsmul`)＋ `Pf` の `k·(α/k) = α` を当てるだけ。 -/
theorem otri_pf_pow (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A)
    (hfn : IsFrobeniusNormalized P A) (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    (ζ : ℕ+ →* End A) (hdeg : ∀ n : ℕ+, P.degFr ((ζ n : End A) : A ⟶ A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P ((ζ n : End A) : A ⟶ A))
    (f : OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :
    ∃ (k : ℕ+) (α : OTri P A),
      ((f : End (⟨A, 1⟩ : PfRootObj P F)) ^ ((k : ℕ+) : ℕ))
        = otriPfMk (F := F) hiso (hprop 1).2 ((α : End A)) := by
  letI : AddCommMonoid (Additive (OTri P A)) := otriAddCommMonoid hfn
  obtain ⟨y, hy⟩ := otriPfMap_surjective (F := F) hiso hA hfn hfn' ζ hdeg hprop f
  induction y using Pf.inductionOn with | _ α k =>
  refine ⟨k, Additive.toMul α, ?_⟩
  have hE := map_nsmul (otriPfEquiv (F := F) hiso hA hfn hfn' ζ hdeg hprop)
    ((k : ℕ+) : ℕ) (Pf.mk α k)
  rw [Pf.nsmul_mk_self] at hE
  have h2 := congrArg (fun t : Additive (OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) =>
    ((Additive.toMul t : OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :
      End (⟨A, 1⟩ : PfRootObj P F))) hE
  rw [← hy]
  exact h2.symm

/-! ## ★5. 根つきと根なしの橋 —— `End ⟨A,1⟩ ≃* Hom^pf(A,A)`

★★`Proposition 5.5, (ii)` の `biratPfHom` の定義域は**根なしの** `Hom^pf(A,A)`
(`PfCat` の `End`)であって、`HomRoot ⟨A,1⟩ ⟨A,1⟩` ではない。
★★★在庫の `endRootMulEquiv`(`rootSelfIso` による `End ⟨A,n⟩ ≃* End_{𝒞^pf} A^{(n·n)}`)を
`n = 1` で読み、`rtExt A 1` が**同型**であること(`isIso_rtExt_one`)による共役を継ぐと橋になる。

★実務メモ: `End (A : PfCat P F)` は圏インスタンスを取り違える(lean-idioms 6)ので、
**別名の `def` `pfObjOf` を通して型を運ぶ**(在庫の `rtObjPf` と同じ手)。 -/

/-- ★`𝒞` の対象を `PfCat` の対象として見る(`End` の解決先を分けるための別名)。 -/
def pfObjOf (P : PreFrobenioid C Φ) (F : FrobenioidCore P) (A : C) : PfCat P F := A

/-- ★★`rtExt A 1` は同型なので、その共役が `End` の同型を与える。 -/
noncomputable def endPfCatRtOne (A : C) :
    End (pfObjOf P F A) ≃* End (pfObjOf P F (rtObj P F A 1)) :=
  Iso.conj ((toPfCat P F).mapIso
    (@asIso _ _ _ _ (rtExt P F A 1) (isIso_rtExt_one P F A)))

/-- ★★★★**根つきと根なしの橋** —— `End ⟨A,1⟩ ≃* Hom^pf(A,A)`。 -/
noncomputable def endRootOneEquiv (A : C) :
    End (⟨A, 1⟩ : PfRootObj P F) ≃* End (pfObjOf P F A) :=
  (endRootMulEquiv (F := F) A 1).trans (endPfCatRtOne (F := F) A).symm

/-- ★★★★★**根なしでも「正の冪で `𝒞` へ戻る」**(`otri_pf_pow` を橋で移したもの)。

★★これが `Proposition 5.5, (ii)` の `biratPfHom` の定義域で使える形である。 -/
theorem hom_pf_pow (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A)
    (hfn : IsFrobeniusNormalized P A) (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    (ζ : ℕ+ →* End A) (hdeg : ∀ n : ℕ+, P.degFr ((ζ n : End A) : A ⟶ A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P ((ζ n : End A) : A ⟶ A))
    (x : End (pfObjOf P F A))
    (hx : (endRootOneEquiv (F := F) A).symm x
      ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :
    ∃ (k : ℕ+) (α : OTri P A),
      x ^ ((k : ℕ+) : ℕ)
        = endRootOneEquiv (F := F) A (otriPfMk (F := F) hiso (hprop 1).2 ((α : End A))) := by
  obtain ⟨k, α, h⟩ := otri_pf_pow hiso hA hfn hfn' ζ hdeg hprop
    ⟨(endRootOneEquiv (F := F) A).symm x, hx⟩
  refine ⟨k, α, ?_⟩
  have hx' : x = endRootOneEquiv (F := F) A ((endRootOneEquiv (F := F) A).symm x) :=
    ((endRootOneEquiv (F := F) A).apply_symm_apply x).symm
  rw [hx', ← map_pow]
  exact congrArg (endRootOneEquiv (F := F) A) h

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Proposition 5.5, (i)` の零因子の側
(`𝒞^pf` の `𝒪^▷` の零因子は `Φ` の元を `ℕ≥1` で割ったもの)。
★これが `Proposition 5.3` の `ℚ·Φ^birat` の `ℚ` の出どころである。 -/
def rootDiv_otri_image.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — Div(𝒪^▷(A^pf)) は Div(𝒪^▷(A)) の ℚ≥0-スパン",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
