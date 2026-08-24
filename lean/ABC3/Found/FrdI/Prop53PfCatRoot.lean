/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Frob

/-!
# [FrdI] `PfCat`(根 1 の部分)から `PfRootObj`(原文の `𝒞^pf`)への関手

原文 (FrdI p.57):
> Also, we obtain a natural functor C → Cpf [by mapping

★★我々の木では `𝒞^pf` を 2 段で作ってある:

* `PfCat P F` —— **根 1 の部分**(対象は `𝒞` の対象、射は `Hom^pf`)
* `PfRootObj P F` —— **原文の `𝒞^pf`**(対象は対 `(A,n)`)

★`𝒞 → PfCat`(`toPfCat`)と `𝒞 → PfRootObj`(`toPfRoot`)は在庫にあるが、
**`PfCat → PfRootObj`** が無かった。★`Proposition 5.3` の図式で
`𝒞^pf` から `𝒞^rlf` へ渡るときに要る(model Frobenioid との圏同値は
`PfRootObj` の側で立ててあるため)。

## ★中身

`rtExt A 1 : A ⟶ A^{(1)}` は**同型**(`isIso_rtExt_one`)なので、
`Hom^pf(A,B)` を `Hom^pf(A^{(1)}, B^{(1)}) = Hom_{𝒞^pf}((A,1),(B,1))` へ
**共役で移す**だけである。

★★関手性は「根 1 どうしでは `compRoot` は `compPf` そのもの」
(`compRoot_root_one`、`compRoot_eq_lift` ＋ `rtRootIso_*_eq_self`)に落ちる。

-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section PfCatRoot

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-- ★★**根 1 の対象どうしでは `compRoot` は `compPf` そのもの**。

★`compRoot_eq_lift` を `c = PA = PB = PE = 1` で当て、
`rtRootIso_hom_eq_self` / `rtRootIso_inv_eq_self`(在庫)で 3 つの輸送を消す。 -/
theorem compRoot_root_one {A B E : C}
    (f : HomRoot P F (⟨A, 1⟩ : PfRootObj P F) ⟨B, 1⟩)
    (g : HomRoot P F (⟨B, 1⟩ : PfRootObj P F) ⟨E, 1⟩) :
    compRoot P F f g = compPf P F f g := by
  have hlift := compRoot_eq_lift (P := P) (F := F) f g
    (c := 1) (PA := 1) (PB := 1) (PE := 1) (hcA := rfl) (hcB := rfl) (hcE := rfl)
    (ef := 1) (eg := 1) (er := 1) (hfA := rfl) (hfB := rfl) (hgA := rfl) (hgE := rfl)
    (hrA := rfl) (hrE := rfl)
  rw [hlift, rtRootIso_inv_eq_self, rtRootIso_inv_eq_self, rtRootIso_hom_eq_self]

/-- ★`PfCat` の対象を `𝒞` の対象として取り出す(型の同義語をほどく)。 -/
abbrev pfObjDown (A : PfCat P F) : C := A

variable {P F} in
/-- ★`PfCat` の合成は `compPf`(`rw` のための橋)。 -/
theorem pfCat_comp_eq {A B E : PfCat P F} (f : A ⟶ B) (g : B ⟶ E) :
    f ≫ g = compPf P F f g := rfl

variable {P F} in
/-- ★`toPfCat` の射の部分は `toHomPf`(`rw` のための橋)。 -/
theorem toPfCat_map_eq {A B : C} (ψ : A ⟶ B) :
    (toPfCat P F).map ψ = toHomPf (F := F) ψ := rfl

/-- ★共役の合成(任意の圏での代数)—— `j ≫ i' = 𝟙` なら真ん中が消える。 -/
theorem conj_comp_aux {X : Type u2} [Category.{v2} X] {a b c a' b' c' : X}
    (i : a' ⟶ a) (f : a ⟶ b) (j : b ⟶ b') (i' : b' ⟶ b) (g : b ⟶ c) (k : c ⟶ c')
    (h : j ≫ i' = 𝟙 b) :
    i ≫ (f ≫ g) ≫ k = (i ≫ f ≫ j) ≫ (i' ≫ g ≫ k) := by
  simp only [Category.assoc]
  rw [← Category.assoc j, h, Category.id_comp]

/-- ★★★★**`PfCat P F ⥤ PfRootObj P F`** —— 根 1 の部分を原文の `𝒞^pf` へ入れる。

★`rtExt A 1` が同型なので、共役で移すだけである。 -/
noncomputable def pfCatToRoot : PfCat P F ⥤ PfRootObj P F where
  obj A := ⟨(pfObjDown P F A), 1⟩
  map {A B} f :=
    haveI := isIso_rtExt_one P F (pfObjDown P F A)
    (show HomRoot P F (⟨(pfObjDown P F A), 1⟩ : PfRootObj P F) ⟨(pfObjDown P F B), 1⟩ from
      (toPfCat P F).map (inv (rtExt P F (pfObjDown P F A) 1)) ≫ f
        ≫ (toPfCat P F).map (rtExt P F (pfObjDown P F B) 1))
  map_id A := by
    haveI := isIso_rtExt_one P F (pfObjDown P F A)
    show (toPfCat P F).map (inv (rtExt P F (pfObjDown P F A) 1))
        ≫ 𝟙 _ ≫ (toPfCat P F).map (rtExt P F (pfObjDown P F A) 1)
      = idRoot P F (⟨(pfObjDown P F A), 1⟩ : PfRootObj P F)
    rw [Category.id_comp, ← (toPfCat P F).map_comp, IsIso.inv_hom_id, (toPfCat P F).map_id]
    rfl
  map_comp {A B E} f g := by
    haveI := isIso_rtExt_one P F (pfObjDown P F A)
    haveI := isIso_rtExt_one P F (pfObjDown P F B)
    haveI := isIso_rtExt_one P F (pfObjDown P F E)
    refine Eq.trans ?_ (compRoot_root_one P F
      ((toPfCat P F).map (inv (rtExt P F (pfObjDown P F A) 1)) ≫ f
        ≫ (toPfCat P F).map (rtExt P F (pfObjDown P F B) 1))
      ((toPfCat P F).map (inv (rtExt P F (pfObjDown P F B) 1)) ≫ g
        ≫ (toPfCat P F).map (rtExt P F (pfObjDown P F E) 1))).symm
    have hcancel : (toPfCat P F).map (rtExt P F (pfObjDown P F B) 1)
        ≫ (toPfCat P F).map (inv (rtExt P F (pfObjDown P F B) 1))
      = 𝟙 B := by
      rw [← (toPfCat P F).map_comp, IsIso.hom_inv_id, (toPfCat P F).map_id]
      rfl
    exact conj_comp_aux _ f _ _ g _ hcancel

variable {P F} in
/-- ★`𝒞 → PfCat → PfRootObj` は `𝒞 → PfRootObj` に等しい(射の側)。 -/
theorem toPfCat_comp_pfCatToRoot_map {A B : C} (φ : A ⟶ B) :
    (toPfCat P F ⋙ pfCatToRoot P F).map φ = (toPfRoot P F).map φ := by
  haveI := isIso_rtExt_one P F A
  haveI := isIso_rtExt_one P F B
  show (toPfCat P F).map (inv (rtExt P F A 1)) ≫ (toPfCat P F).map φ
      ≫ (toPfCat P F).map (rtExt P F B 1)
    = toRootHom (F := F) φ
  rw [← (toPfCat P F).map_comp, ← (toPfCat P F).map_comp]
  rfl

/-! ## ★2. 残り(記録)

★★**`𝒞^pf`(根 1 の部分 `PfCat`)が isotropic 型であること**は未着手である。
在庫の `pfRoot_isOfIsotropicType`(`PfRootObj` の言葉、`Proposition 3.2, (iii)`)を
`pfCatToRoot` で降ろせばよい —— `pfDiv`/`pfDeg`/`pfBase` が共役で不変であること
(`rtExt A 1` は Frobenius 型の同型なので `Div = 0`、`deg = 1`、`Base` 同型)と、
`compRoot_root_one` で同型性を `PfCat` へ戻すことの 2 段。

★★**2026-08-24 の試作で分かったこと(次に着手する人へ)**:

1. 詰まりの主因は **`pfDown` を使わずに自前の `pfObjDown` を書いたこと**だった。
   `pfDiv` などの暗黙引数は在庫の `pfDown` で書かれているので、
   **在庫の `pfDown` に揃える**と `pfDiv_comp` / `pfDeg_comp` /
   `rootDiv` / `rootBase` / `compRoot_root_one` はすべて当たる。
2. `≫` を `compPf` に開く橋 `pfCat_comp_eq` と、
   `(toPfCat P F).map ψ = toHomPf ψ` の橋 `toPfCat_map_eq` を用意し、
   **目標の側で `rw` せず、自分で述べた `have` の側で `rw` する**
   (`≫` を開いた目標は `instances` 透明度で型が付かなくなる)。
3. それでも残ったのは **インスタンス合成 3 件**:
   `IsIso (P.Base (inv (rtExt X 1)) ≫ pfBase φ ≫ P.Base (rtExt Dd 1))`、
   `φ = φ ≫ j ≫ inv j` の後始末、最後の `IsIso (inv i ≫ … ≫ inv j)`。
   ★どれも「暗黙引数の書かれ方が違うので既存の `IsIso` 仮定が拾われない」型。

★`Proposition 5.3` の図式ではこれを `hisoPf` として仮引数で受けている。
-/

end PfCatRoot

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Definition 3.1, (iii)` の「根 1 の部分から `𝒞^pf` へ」。 -/
def pfCatToRoot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 57,
    item := "Definition 3.1, (iii) — 根 1 の部分 PfCat から 𝒞^pf(PfRootObj)への関手",
    sectionId := "frdi-def-3-1" }

end ABC3.Found.FrdI
