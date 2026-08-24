/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55Pf
import ABC3.Found.FrdI.Prop32Equiv
import ABC3.Found.FrdI.Prop53UntrPf
import ABC3.Found.FrdI.Def24PfT
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

/-- ★★★★★**同上を「自然な関手の像」の形で** —— 添字 1 の `otriPfMk` は
`toRootHom`(自然な関手 `𝒞 → 𝒞^pf` が定める写像)そのもの(`otriPfMk_one`)。 -/
theorem otri_pf_pow_toRootHom (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A)
    (hfn : IsFrobeniusNormalized P A) (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    (ζ : ℕ+ →* End A) (hdeg : ∀ n : ℕ+, P.degFr ((ζ n : End A) : A ⟶ A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P ((ζ n : End A) : A ⟶ A))
    (f : OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :
    ∃ (k : ℕ+) (α : OTri P A),
      ((f : End (⟨A, 1⟩ : PfRootObj P F)) ^ ((k : ℕ+) : ℕ))
        = toRootHom (F := F) ((α : End A) : A ⟶ A) := by
  obtain ⟨k, α, h⟩ := otri_pf_pow hiso hA hfn hfn' ζ hdeg hprop f
  exact ⟨k, α, h.trans (otriPfMk_one hiso (hprop 1).2
    (congrArg (fun t : End A => (t : A ⟶ A)) (map_one ζ)) ((α : End A)))⟩

/-- ★★★★★同上、根なし版。 -/
theorem hom_pf_pow_toRootHom (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A)
    (hfn : IsFrobeniusNormalized P A) (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    (ζ : ℕ+ →* End A) (hdeg : ∀ n : ℕ+, P.degFr ((ζ n : End A) : A ⟶ A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P ((ζ n : End A) : A ⟶ A))
    (x : End (pfObjOf P F A))
    (hx : (endRootOneEquiv (F := F) A).symm x
      ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :
    ∃ (k : ℕ+) (α : OTri P A),
      x ^ ((k : ℕ+) : ℕ)
        = endRootOneEquiv (F := F) A (toRootHom (F := F) ((α : End A) : A ⟶ A)) := by
  obtain ⟨k, α, h⟩ := otri_pf_pow_toRootHom hiso hA hfn hfn' ζ hdeg hprop
    ⟨(endRootOneEquiv (F := F) A).symm x, hx⟩
  refine ⟨k, α, ?_⟩
  have hx' : x = endRootOneEquiv (F := F) A ((endRootOneEquiv (F := F) A).symm x) :=
    ((endRootOneEquiv (F := F) A).apply_symm_apply x).symm
  rw [hx', ← map_pow]
  exact congrArg (endRootOneEquiv (F := F) A) h

/-- ★★`n = 1` では `Θ`(`rootSelfIso`)は恒等 —— `rtRootIso` の `e = 1` の場合。 -/
theorem rootSelfIso_one_inv (A : C) (z : HomPf P F (rtObj P F A 1) (rtObj P F A 1)) :
    (rootSelfIso (F := F) A 1).inv z = z :=
  rtRootIso_inv_eq_self P F A A (show (1 : ℕ+) = 1 * 1 from rfl)
    (show (1 : ℕ+) = 1 * 1 from rfl) z

/-- ★★★★**橋は自然な関手と両立する** —— `endRootOneEquiv (toRootHom α) = toHomPf α`。

★`Θ` は `n = 1` では恒等、残る共役は `rtExt A 1` とその逆で打ち消し合う。 -/
theorem endRootOneEquiv_toRootHom (A : C) (α : A ⟶ A) :
    endRootOneEquiv (F := F) A (toRootHom (F := F) α) = toHomPf (F := F) α := by
  haveI := isIso_rtExt_one P F A
  show (endPfCatRtOne (F := F) A).symm ((rootSelfIso (F := F) A 1).inv (toRootHom (F := F) α))
      = toHomPf (F := F) α
  rw [rootSelfIso_one_inv]
  show compPf P F (toHomPf (F := F) (rtExt P F A 1))
      (compPf P F (toRootHom (F := F) α) (toHomPf (F := F) (inv (rtExt P F A 1))))
    = toHomPf (F := F) α
  show compPf P F (toHomPf (F := F) (rtExt P F A 1))
      (compPf P F (toHomPf (F := F) (inv (rtExt P F A 1) ≫ α ≫ rtExt P F A 1))
        (toHomPf (F := F) (inv (rtExt P F A 1))))
    = toHomPf (F := F) α
  rw [← toHomPf_comp, ← toHomPf_comp]
  refine congrArg (toHomPf (F := F)) ?_
  simp

/-! ## ★6. 逆に —— `𝒪^▷(A^pf)` は**可除**である

★★★これが `(Φ^pf)^birat = ℚ·Φ^birat` の**等号**に要る側である。
`𝒪^▷(A^pf) ≅ 𝒪^▷(A)^pf` の下で `k · (α/(m·k)) = α/m` を読み替えるだけ。 -/

/-- ★`Pf` の中の可除性 —— `k · (a/(m·k)) = a/m`。 -/
theorem Pf.nsmul_mk_mul {M : Type w} [AddCommMonoid M] (a : M) (m k : ℕ+) :
    ((k : ℕ+) : ℕ) • Pf.mk a (m * k) = Pf.mk a m := by
  rw [Pf.nsmul_mk]
  refine (Pf.mk_nsmul_eq_of_cross' (N := ((k : ℕ+) : ℕ)) (D := m * k) (P := 1) (Q := m) a ?_).trans
    (by rw [one_smul])
  push_cast
  ring

/-- ★★★★★**`𝒪^▷(A^pf)` は可除** —— どの元にも `k` 乗根が `𝒪^▷(A^pf)` の中にある。 -/
theorem otri_pf_root (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A)
    (hfn : IsFrobeniusNormalized P A) (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    (ζ : ℕ+ →* End A) (hdeg : ∀ n : ℕ+, P.degFr ((ζ n : End A) : A ⟶ A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P ((ζ n : End A) : A ⟶ A))
    (f : OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) (k : ℕ+) :
    ∃ w : OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F),
      ((w : End (⟨A, 1⟩ : PfRootObj P F)) ^ ((k : ℕ+) : ℕ))
        = (f : End (⟨A, 1⟩ : PfRootObj P F)) := by
  letI : AddCommMonoid (Additive (OTri P A)) := otriAddCommMonoid hfn
  obtain ⟨y, hy⟩ := otriPfMap_surjective (F := F) hiso hA hfn hfn' ζ hdeg hprop f
  induction y using Pf.inductionOn with | _ α m =>
  refine ⟨otriPfMap (F := F) hiso hfn ζ hdeg hprop (Pf.mk α (m * k)), ?_⟩
  have hE := map_nsmul (otriPfEquiv (F := F) hiso hA hfn hfn' ζ hdeg hprop)
    ((k : ℕ+) : ℕ) (Pf.mk α (m * k))
  rw [Pf.nsmul_mk_mul] at hE
  have h2 := congrArg (fun t : Additive (OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) =>
    ((Additive.toMul t : OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :
      End (⟨A, 1⟩ : PfRootObj P F))) hE
  rw [← hy]
  exact h2.symm

/-- ★★★★★同上、根なし版(`biratPfHom` の定義域で使える形)。 -/
theorem hom_pf_root (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A)
    (hfn : IsFrobeniusNormalized P A) (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    (ζ : ℕ+ →* End A) (hdeg : ∀ n : ℕ+, P.degFr ((ζ n : End A) : A ⟶ A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P ((ζ n : End A) : A ⟶ A))
    (x : End (pfObjOf P F A))
    (hx : (endRootOneEquiv (F := F) A).symm x
      ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) (k : ℕ+) :
    ∃ x' : End (pfObjOf P F A),
      (endRootOneEquiv (F := F) A).symm x' ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)
        ∧ x' ^ ((k : ℕ+) : ℕ) = x := by
  obtain ⟨w, hw⟩ := otri_pf_root hiso hA hfn hfn' ζ hdeg hprop
    ⟨(endRootOneEquiv (F := F) A).symm x, hx⟩ k
  refine ⟨endRootOneEquiv (F := F) A ((w : End (⟨A, 1⟩ : PfRootObj P F))), ?_, ?_⟩
  · rw [(endRootOneEquiv (F := F) A).symm_apply_apply]
    exact w.2
  · rw [← map_pow, hw]
    exact (endRootOneEquiv (F := F) A).apply_symm_apply x

/-! ## ★7. `Gp (Pf M)` は torsion-free

★★★これが「飽和 = それ自身」を言うのに要る最後の道具である。
`Pf M` では `n •` が**全単射**(単射は在庫 `Pf.nsmul_injective`、
全射は `Pf.nsmul_mk_mul`)なので、群化しても捩れは出ない。 -/

/-- ★★`Pf M` では `n •` は全射(`a/m = n · (a/(m·n))`)。 -/
theorem Pf.nsmul_surjective {M : Type w} [AddCommMonoid M] (k : ℕ+) :
    Function.Surjective (fun x : Pf M => ((k : ℕ+) : ℕ) • x) := by
  intro x
  induction x using Pf.inductionOn with | _ m a =>
  exact ⟨Pf.mk m (a * k), Pf.nsmul_mk_mul m a k⟩

/-- ★★★★**`Gp (Pf M)` は torsion-free**(仮定なし)。 -/
theorem gp_pf_isTorsionFreeNaive {M : Type w} [AddCommMonoid M] :
    IsTorsionFreeNaive (Gp (Pf M)) := by
  intro x n hn hx
  obtain ⟨a, b, rfl⟩ : ∃ (a b : Pf M), x = toGp (Pf M) a - toGp (Pf M) b := by
    induction x using AddLocalization.induction_on with
    | H y => exact ⟨y.1, (y.2 : Pf M), (eq_sub_of_add_eq (mk_add_toGp (Pf M) y.1 y.2)).symm ▸ rfl⟩
  rw [smul_sub, ← toGp_nsmul, ← toGp_nsmul, sub_eq_zero] at hx
  obtain ⟨c, hc⟩ := toGp_eq_iff.mp hx
  obtain ⟨c', rfl⟩ := Pf.nsmul_surjective (M := M) ⟨n, hn⟩ c
  have hc2 : n • (a + c') = n • (b + c') := by
    rw [smul_add, smul_add]
    exact hc
  have hab := Pf.nsmul_injective (M := M) ⟨n, hn⟩ hc2
  rw [sub_eq_zero]
  exact toGp_eq_iff.mpr ⟨c', hab⟩

/-- ★★**同型による共役は `𝒪^▷` を保つ**(底恒等性も次数 1 も共役で不変)。 -/
theorem mem_otri_conj {Cc : Type u2} [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} D}
    (P₀ : PreFrobenioid Cc Φ₀) {X Y : Cc} (u : X ≅ Y) (f : End X)
    (hf : f ∈ OTri P₀ X) : (u.conj f) ∈ OTri P₀ Y := by
  obtain ⟨hb, hl⟩ := hf
  haveI : IsIso u.hom := u.isIso_hom
  haveI : IsIso u.inv := u.isIso_inv
  refine ⟨?_, ?_⟩
  · show P₀.Base (u.inv ≫ ((f : End X) : X ⟶ X) ≫ u.hom) = P₀.Base (𝟙 Y)
    rw [P₀.Base_comp, P₀.Base_comp, show P₀.Base ((f : End X) : X ⟶ X) = P₀.Base (𝟙 X) from hb,
      P₀.Base_id, Category.id_comp, ← P₀.Base_comp, u.inv_hom_id]
  · show P₀.degFr (u.inv ≫ ((f : End X) : X ⟶ X) ≫ u.hom) = 1
    rw [P₀.degFr_comp, P₀.degFr_comp, degFr_of_isIso P₀ u.hom, degFr_of_isIso P₀ u.inv,
      show P₀.degFr ((f : End X) : X ⟶ X) = 1 from hl, one_mul, one_mul]

/-- ★★★**根なしの `𝒪^▷` は根つき化身の `𝒪^▷` に移る**(橋 `endRootOneEquiv` で)。 -/
theorem mem_otri_endRootOneEquiv_symm (A : C) (x : End (pfObjOf P F A))
    (hx : x ∈ OTri (pfPre P F) (pfObjOf P F A)) :
    (endRootOneEquiv (F := F) A).symm x
      ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F) := by
  refine (mem_otri_rootSelfIso_inv A 1 _).mpr ?_
  have h : endRootMulEquiv (F := F) A 1 ((endRootOneEquiv (F := F) A).symm x)
      = endPfCatRtOne (F := F) A x :=
    (endRootMulEquiv (F := F) A 1).apply_symm_apply _
  show endRootMulEquiv (F := F) A 1 ((endRootOneEquiv (F := F) A).symm x)
      ∈ OTri (pfPre P F) (rtObjPf P F A 1)
  rw [h]
  exact mem_otri_conj (pfPre P F)
    ((toPfCat P F).mapIso (@asIso _ _ _ _ (rtExt P F A 1) (isIso_rtExt_one P F A))) x hx

/-- ★★同型による共役は `𝒪^▷` を**反射**もする(`mem_otri_conj` を逆向きの同型に当てる)。 -/
theorem mem_otri_endRootOneEquiv_symm_of (A : C) (x : End (pfObjOf P F A))
    (hx : (endRootOneEquiv (F := F) A).symm x
      ∈ OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F)) :
    x ∈ OTri (pfPre P F) (pfObjOf P F A) := by
  have h : endRootMulEquiv (F := F) A 1 ((endRootOneEquiv (F := F) A).symm x)
      = endPfCatRtOne (F := F) A x :=
    (endRootMulEquiv (F := F) A 1).apply_symm_apply _
  have h1 := (mem_otri_rootSelfIso_inv A 1 _).mp hx
  have h2 : endPfCatRtOne (F := F) A x ∈ OTri (pfPre P F) (rtObjPf P F A 1) := by
    rw [← h]
    exact h1
  have h3 := mem_otri_conj (pfPre P F)
    ((toPfCat P F).mapIso (@asIso _ _ _ _ (rtExt P F A 1) (isIso_rtExt_one P F A))).symm
    (endPfCatRtOne (F := F) A x) h2
  have h4 : (((toPfCat P F).mapIso
        (@asIso _ _ _ _ (rtExt P F A 1) (isIso_rtExt_one P F A))).symm.conj
      (endPfCatRtOne (F := F) A x)) = x :=
    (endPfCatRtOne (F := F) A).symm_apply_apply x
  rw [h4] at h3
  exact h3

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Proposition 5.5, (i)` の零因子の側
(`𝒞^pf` の `𝒪^▷` の零因子は `Φ` の元を `ℕ≥1` で割ったもの)。
★これが `Proposition 5.3` の `ℚ·Φ^birat` の `ℚ` の出どころである。 -/
def rootDiv_otri_image.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — Div(𝒪^▷(A^pf)) は Div(𝒪^▷(A)) の ℚ≥0-スパン",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
