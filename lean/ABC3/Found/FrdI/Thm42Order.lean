/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42DivMap
import ABC3.Found.FrdI.Thm42PsiPrime

/-!
# [FrdI] Theorem 4.2, (ii) の `Order` と (iii) の単系同型

原文 (FrdI p.77):
> for the respective full subcategories and restricted equivalences of categories

原文 (FrdI p.78):
> then the last two equivalences of categories of (ii) arise from isomorphisms of

## ★★★★★測って分かった —— **単系同型が先、圏同値は後**

原文の順序は「圏同値 (ii) → (Div-Frobenius-trivial なら) 単系同型 (iii)」だが、
★我々は `divEquiv`(`Thm42DivMap.lean`、逸脱 (B) `hdivS` を使う)で
**先に単系同型 `Φ(A) ≅ Φ₂(ΨA)` を持っている**。

★★したがって順序が逆になる:
1. `divEquiv` を `M_p` に**制限**する(`mpAddEquiv`)——★素点の対応は
   `primeEquiv (divEquiv …)` であり、これが **`psiPrime` と一致する**
   (`psiPrime_eq_primeEquiv`)。
2. `Order(-)` は前順序圏なので、単系同型はそのまま**圏同値**を与える
   (`orderCatEquivalence`)。反対圏版は `Equivalence.op`。

★★★**逸脱ではない** —— 結論は原文と同じで、依存の向きだけが違う。
`hdivS` は `prop_4_4.src` で開示済みの逸脱 (B) である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★1. `Order(-)` は単系同型で移る -/

section OrderTransport

variable {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]

/-- ★加法準同型が定める `Order(M) ⥤ Order(N)`。 -/
def orderCatFunctor (f : M →+ N) : OrderCat M ⥤ OrderCat N where
  obj a := toOrderCat (f (a : M))
  map {a b} h := homOfLE (MLe.map f (leOfHom h))

/-- ★★**単系同型は `Order(-)` の圏同値を与える** ——
`Order(-)` は thin な圏なので、対象の対応さえ付けば自然性は自動。 -/
def orderCatEquivalence (e : M ≃+ N) : OrderCat M ≌ OrderCat N where
  functor := orderCatFunctor e.toAddMonoidHom
  inverse := orderCatFunctor e.symm.toAddMonoidHom
  unitIso := NatIso.ofComponents
    (fun a => eqToIso (show a = toOrderCat (e.symm (e (a : M))) by
      rw [e.symm_apply_apply]; rfl))
    (fun _ => Subsingleton.elim _ _)
  counitIso := NatIso.ofComponents
    (fun b => eqToIso (show toOrderCat (e (e.symm (b : N))) = b by
      rw [e.apply_symm_apply]; rfl))
    (fun _ => Subsingleton.elim _ _)
  functor_unitIso_comp _ := Subsingleton.elim _ _

def orderCatEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 12,
    item := "§0 Monoids — Order(M) は単系同型で移る",
    sectionId := "frdi-s0-prime" }

/-! ## ★2. `M_p` の対応 -/

/-- ★`M_p` の元は `e` で `N_{e(p)}` へ移る。 -/
theorem Mp_map_mem (e : M ≃+ N) (p : Prime M) {a : M} (h : a ∈ Mp M p) :
    e a ∈ Mp N (primeEquiv e p) := by
  refine AddSubmonoid.closure_induction ?_ ?_ ?_ h
  · intro x hx
    exact mem_Mp_of_mem_primeCarrier ((mem_primeCarrier_map e p x).mp hx)
  · rw [map_zero]
    exact zero_mem _
  · intro x y _ _ hx hy
    rw [map_add]
    exact add_mem hx hy

/-- ★★`a ∈ M_p ⟺ e a ∈ N_{e(p)}`。 -/
theorem mem_Mp_map (e : M ≃+ N) (p : Prime M) (a : M) :
    a ∈ Mp M p ↔ e a ∈ Mp N (primeEquiv e p) := by
  refine ⟨Mp_map_mem e p, fun h => ?_⟩
  have h2 := Mp_map_mem e.symm (primeEquiv e p) h
  rw [e.symm_apply_apply] at h2
  have h3 : primeEquiv e.symm (primeEquiv e p) = p := (primeEquiv e).left_inv p
  rwa [h3] at h2

/-- ★★★**単系同型は `M_p ≅ N_{e(p)}` を与える**。 -/
def mpAddEquiv (e : M ≃+ N) (p : Prime M) : Mp M p ≃+ Mp N (primeEquiv e p) where
  toFun x := ⟨e x.1, (mem_Mp_map e p x.1).mp x.2⟩
  invFun y := ⟨e.symm y.1, (mem_Mp_map e p (e.symm y.1)).mpr
    (by rw [e.apply_symm_apply]; exact y.2)⟩
  left_inv x := Subtype.ext (e.symm_apply_apply x.1)
  right_inv y := Subtype.ext (e.apply_symm_apply y.1)
  map_add' x y := Subtype.ext (map_add e x.1 y.1)

def mpAddEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 12,
    item := "§0 Monoids — M_p は単系同型で移る",
    sectionId := "frdi-s0-prime" }

/-- ★★★★**`Order(M_p) ≌ Order(N_{e(p)})`**。 -/
def orderMpEquivalence (e : M ≃+ N) (p : Prime M) :
    OrderCat (Mp M p) ≌ OrderCat (Mp N (primeEquiv e p)) :=
  orderCatEquivalence (mpAddEquiv e p)

/-- ★★★★**反対圏版**。 -/
def orderMpEquivalenceOp (e : M ≃+ N) (p : Prime M) :
    (OrderCat (Mp M p))ᵒᵖ ≌ (OrderCat (Mp N (primeEquiv e p)))ᵒᵖ :=
  (orderMpEquivalence e p).op

end OrderTransport

/-! ## ★3. `Ψ_Prime` と `divEquiv` の整合 -/

section PsiOrder

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

variable (G : Frobenioid P) (G₂ : Frobenioid P₂) (Ψ : C ≌ C₂)

variable (P P₂) in
/-- `Ψ` が `𝒪^▷` を保つ条。 -/
abbrev OTriFwd : Prop := ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
  ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
    ∈ OTri P₂ (Ψ.functor.obj Z)

variable (P P₂) in
/-- `Ψ` が `𝒪^▷` を反射する条。 -/
abbrev OTriBwd : Prop := ∀ (Z : C) (δ : End Z),
  ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
    ∈ OTri P₂ (Ψ.functor.obj Z) → δ ∈ OTri P Z

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`Ψ_Prime` は `divEquiv` が誘導する素点の対応に他ならない**。

★これが (ii) の `Order` と (iii) の単系同型を結ぶ蝶番である。 -/
theorem psiPrime_eq_primeEquiv (ctx : PrimeCtx P P₂ G G₂ Ψ)
    (hOTri : OTriFwd P P₂ Ψ) (hOTri' : OTriBwd P P₂ Ψ)
    (A : C) (p : Prime (Φ.val (P.toElem.obj A).base)) :
    psiPrime ctx A p
      = primeEquiv (divEquiv Ψ G G₂ ctx.iso ctx.iso₂ ctx.divS ctx.divS₂ hOTri hOTri' A) p := by
  -- `u` は `p` を実現する `𝒪^▷` の自己射
  set u : End A := realizeIn ctx A p with hu
  have hum : u ∈ OTri P A := realizeIn_mem ctx A p
  have hus : IsPreStep P (((u : End A)) : A ⟶ A) := realizeIn_preStep ctx A p
  have hupv : preStepVal P (((u : End A)) : A ⟶ A) hus = P.Div ((((u : End A)) : A ⟶ A)) :=
    preStepVal_of_otri u hum hus
  -- `p = toPrime (Div u)`
  have hp0 : toPrime _ (preStepVal P (((u : End A)) : A ⟶ A) hus) (realizeIn_primary ctx A p)
      = p := realizeIn_toPrime ctx A p
  have hprimU : IsPrimaryElt (P.Div ((((u : End A)) : A ⟶ A))) := by
    rw [← hupv]; exact realizeIn_primary ctx A p
  have hp : toPrime _ (P.Div ((((u : End A)) : A ⟶ A))) hprimU = p := by
    rw [← hp0]; exact toPrime_congr hupv.symm
  -- `Ψ u` の側
  have hvm : ((Ψ.functor.map ((((u : End A)) : A ⟶ A))) : End (Ψ.functor.obj A))
      ∈ OTri P₂ (Ψ.functor.obj A) := hOTri A u hum
  have hvs : IsPreStep P₂ (Ψ.functor.map ((((u : End A)) : A ⟶ A))) := ctx.PS _ hus
  have hvpv : preStepVal P₂ (Ψ.functor.map ((((u : End A)) : A ⟶ A))) hvs
      = P₂.Div (Ψ.functor.map ((((u : End A)) : A ⟶ A))) :=
    preStepVal_of_otri _ hvm hvs
  -- ★`divEquiv (Div u) = Div (Ψ u)`
  have hde : divEquiv Ψ G G₂ ctx.iso ctx.iso₂ ctx.divS ctx.divS₂ hOTri hOTri' A
      (P.Div ((((u : End A)) : A ⟶ A)))
      = P₂.Div (Ψ.functor.map ((((u : End A)) : A ⟶ A))) := by
    show P₂.Div (Ψ.functor.map ((((realizeDiv P ctx.divS A
      (P.Div ((((u : End A)) : A ⟶ A)))) : End A) : A ⟶ A))) = _
    exact div_map_eq_of_div_eq Ψ G ctx.iso
      (realizeDiv_mem ctx.divS A (P.Div ((((u : End A)) : A ⟶ A)))) hum
      (by rw [realizeDiv_div])
  -- 両辺を `toPrime (Div (Ψ u))` に落とす
  have hl : psiPrime ctx A p
      = toPrime _ (P₂.Div (Ψ.functor.map ((((u : End A)) : A ⟶ A))))
          (by rw [← hvpv]; exact map_realizeIn_primary ctx A p) := by
    show toPrime _ (preStepVal P₂ (Ψ.functor.map ((((u : End A)) : A ⟶ A))) _) _ = _
    exact toPrime_congr hvpv
  have hr : primeEquiv (divEquiv Ψ G G₂ ctx.iso ctx.iso₂ ctx.divS ctx.divS₂ hOTri hOTri' A) p
      = toPrime _ (P₂.Div (Ψ.functor.map ((((u : End A)) : A ⟶ A))))
          (by rw [← hde]; exact isPrimaryElt_map _ hprimU) := by
    rw [← hp]
    show toPrime _ (divEquiv Ψ G G₂ ctx.iso ctx.iso₂ ctx.divS ctx.divS₂ hOTri hOTri' A
      (P.Div ((((u : End A)) : A ⟶ A)))) _ = _
    exact toPrime_congr hde
  rw [hl, hr]

def psiPrime_eq_primeEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (iii) — Ψ_Prime は divEquiv が誘導する素点の対応",
    sectionId := "frdi-thm-4-2" }

/-! ## ★4. 結論 -/

/-- ★★★★★★**[FrdI] Theorem 4.2, (iii)** ——
`Ψ` は**単系同型 `Φ(A)_p ≅ Φ₂(ΨA)_{Ψ_Prime p}`** を誘導する。 -/
noncomputable def thm_4_2_iii_monoid (ctx : PrimeCtx P P₂ G G₂ Ψ)
    (hOTri : OTriFwd P P₂ Ψ) (hOTri' : OTriBwd P P₂ Ψ)
    (A : C) (p : Prime (Φ.val (P.toElem.obj A).base)) :
    Mp (Φ.val (P.toElem.obj A).base) p
      ≃+ Mp (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base) (psiPrime ctx A p) := by
  rw [psiPrime_eq_primeEquiv G G₂ Ψ ctx hOTri hOTri' A p]
  exact mpAddEquiv _ p

def thm_4_2_iii_monoid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (iii)", sectionId := "frdi-thm-4-2" }

/-- ★★★★★★**[FrdI] Theorem 4.2, (ii)** ——
`Order(Φ(A)_p) ≌ Order(Φ₂(ΨA)_{Ψ_Prime p})`。 -/
noncomputable def thm_4_2_ii_order (ctx : PrimeCtx P P₂ G G₂ Ψ)
    (hOTri : OTriFwd P P₂ Ψ) (hOTri' : OTriBwd P P₂ Ψ)
    (A : C) (p : Prime (Φ.val (P.toElem.obj A).base)) :
    OrderCat (Mp (Φ.val (P.toElem.obj A).base) p)
      ≌ OrderCat (Mp (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base) (psiPrime ctx A p)) :=
  orderCatEquivalence (thm_4_2_iii_monoid G G₂ Ψ ctx hOTri hOTri' A p)

/-- ★★★★★★**[FrdI] Theorem 4.2, (ii)** —— 反対圏版。 -/
noncomputable def thm_4_2_ii_order_op (ctx : PrimeCtx P P₂ G G₂ Ψ)
    (hOTri : OTriFwd P P₂ Ψ) (hOTri' : OTriBwd P P₂ Ψ)
    (A : C) (p : Prime (Φ.val (P.toElem.obj A).base)) :
    (OrderCat (Mp (Φ.val (P.toElem.obj A).base) p))ᵒᵖ
      ≌ (OrderCat (Mp (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base) (psiPrime ctx A p)))ᵒᵖ :=
  (thm_4_2_ii_order G G₂ Ψ ctx hOTri hOTri' A p).op

def thm_4_2_ii_order.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (ii) — Order(Φ(A)_p) の圏同値",
    sectionId := "frdi-thm-4-2" }

/-! ## ★5. `Ψ_Prime` の一意性 -/

/-- ★★★★★**[FrdI] Theorem 4.2, (ii) の「一意な」** ——
原文の特徴づけ(primary pre-step の対応)を満たす写像は `psiPrime` に限る。

★どの素点も primary pre-step で実現される(`realizeIn`)ので、
特徴づけが写像を完全に決める。 -/
theorem psiPrime_unique (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C)
    (F : Prime (Φ.val (P.toElem.obj A).base) →
      Prime (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base))
    (hF : ∀ {E : C} (ϵ : E ⟶ A) (hsϵ : IsPreStep P ϵ)
      (hpϵ : IsPrimaryElt (preStepVal P ϵ hsϵ))
      (hq : IsPrimaryElt (preStepVal P₂ (Ψ.functor.map ϵ) (ctx.PS ϵ hsϵ))),
      F (toPrime _ (preStepVal P ϵ hsϵ) hpϵ)
        = toPrime _ (preStepVal P₂ (Ψ.functor.map ϵ) (ctx.PS ϵ hsϵ)) hq) :
    F = psiPrime ctx A := by
  funext p
  have hu := realizeIn_toPrime ctx A p
  have h2 := psiPrime_spec ctx (((realizeIn ctx A p) : A ⟶ A)) (realizeIn_preStep ctx A p)
    (realizeIn_primary ctx A p)
  rw [hu] at h2
  have h3 := hF (((realizeIn ctx A p) : A ⟶ A)) (realizeIn_preStep ctx A p)
    (realizeIn_primary ctx A p) (map_realizeIn_primary ctx A p)
  rw [hu] at h3
  rw [h3, h2]

def psiPrime_unique.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (ii) — Ψ_Prime の一意性",
    sectionId := "frdi-thm-4-2" }

/-! ## ★6. `Theorem 4.2` 全体 -/

/-- ★★★★★★★**[FrdI] Theorem 4.2** —— 条なしの locator。

| 主張 | 実装 |
|---|---|
| (i) primary step の保存 | `primaryInIff_all`(`Thm42Prop41v.lean`) |
| (i) Div-identity 自己射の保存 | `isDivIdentity_map`(`Thm42DivId.lean`) |
| (i) (universally) Div-Frobenius-trivial 対象の保存 | `divFrobeniusTrivial_map` |
| (ii) `Ψ_Prime` の存在・全単射・関手性 | `psiPrime` / `psiPrimeEquiv` / `psiPrime_naturality` |
| (ii) 一意性 | `psiPrime_unique` |
| (ii) `Order(Φ(A)_p)` の圏同値(＋反対圏) | `thm_4_2_ii_order` / `thm_4_2_ii_order_op` |
| (iii) 単系同型 | `thm_4_2_iii_monoid` |

★★逸脱 (B)(`hdivS`: `Div : 𝒪^▷(A) ↠ Φ(A)`、`prop_4_4.src` で開示済)を使うので、
(iii) は原文の `A_i` が Div-Frobenius-trivial という条を**要さない**
(原文より強い結論であり、後続の消費に影響しない)。 -/
def thm_4_2.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77, item := "Theorem 4.2",
    sectionId := "frdi-thm-4-2" }

end PsiOrder

end ABC3.Found.FrdI
