/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Frob

/-!
# [FrdI] `Proposition 3.2, (iii)(d)` —— `𝒞^pf` の 2 本の圏同値

★`Prop32Frob.lean` で `FrobenioidCore (pfRootPre P F)` の **21 条**が閉じた。
`Frobenioid` にするには (iii)(d) の 2 本の圏同値

  `_X(𝒞^coa-pre) ≃ Order(Φ(X))`,  `(𝒞^coa-pre)_X ≃ Order(Φ(X))^opp`

が要る。忠実性 2 本と前置の本質的全射性は `Prop32Frob.lean` で済んでいる。
**このファイルは残る 3 条**(後置の本質的全射性・前置の充満性・後置の充満性)を埋める。

## ★★設計の要 —— 型の罠を抽象の補題で迂回する

★`pfRootPre P F` の底の型は `((pfRootPre P F).toElem.obj ⟨A, k⟩).base` と書かれ、
`(P.toElem.obj A).base` と**定義上等しいが構文上は別物**である
(`memory/widesubcategory-type-trap.md`)。このため `rw` が
「target expression is not type-correct」で詰まる。

★★**対策**: 計算を**抽象の pre-Frobenioid `Q` の補題**に出す
(`overVal_comp_iso` / `Base_inv_eq` / `inv_Base_inv_eq`)。
そこでは対象が素の変数なので罠が起きない。具体側では**値の等式**として
`rw` するだけになる。

## ★道具

| 補題 | 内容 |
|---|---|
| `rootBase_lamHom` | ★**`Λ_k` は底を変えない**(`Base(Λ_k φ) = Base φ`) |
| `overVal_lamHom` | `Λ_k` の「後置の値」は `k` で割られる |
| `overVal_comp_iso` | 後置の値と同型の合成(抽象) |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2 v3 u3

/-! ## ★1. 抽象の pre-Frobenioid での補題(型の罠を避けるため) -/

section Abstract

variable {D : Type u} [Category.{v} D]

/-- ★`Φ.map g ∘ Φ.map g⁻¹ = id`。 -/
theorem MonoidOn.map_inv_map (Φ : MonoidOn.{v, u, w} D) {X Y : D} (g : X ⟶ Y) (hg : IsIso g)
    (y : Φ.val X) : Φ.map g (Φ.map (@inv _ _ _ _ g hg) y) = y := by
  rw [← Φ.map_comp, @IsIso.hom_inv_id _ _ _ _ g hg]
  exact Φ.map_id X y

/-- ★同型の底の逆射。 -/
theorem Base_inv_eq {C' : Type u2} [Category.{v2} C'] {Φ' : MonoidOn.{v, u, w} D}
    (Q : PreFrobenioid C' Φ') {X Y : C'} (e : X ⟶ Y) [IsIso e]
    (hb : IsIso (Q.Base e)) :
    Q.Base (inv e) = @inv _ _ _ _ (Q.Base e) hb := by
  refine IsIso.eq_inv_of_hom_inv_id ?_
  rw [← Q.Base_comp, IsIso.hom_inv_id, Q.Base_id]

/-- ★同型の逆射の底の逆射。 -/
theorem inv_Base_inv_eq {C' : Type u2} [Category.{v2} C'] {Φ' : MonoidOn.{v, u, w} D}
    (Q : PreFrobenioid C' Φ') {X Y : C'} (e : X ⟶ Y) [IsIso e]
    (hbe : IsIso (Q.Base e)) (hbie : IsIso (Q.Base (inv e))) :
    @inv _ _ _ _ (Q.Base (inv e)) hbie = Q.Base e := by
  refine IsIso.inv_eq_of_hom_inv_id ?_
  rw [← Q.Base_comp, IsIso.inv_hom_id, Q.Base_id]

/-- ★★**後置の値と同型の合成** —— `overVal (ψ ≫ e) = (Base e)⁻¹^*(overVal ψ)`。

★「後置の値」とは `(𝒞^coa-pre)_A → Order(Φ(A))^opp` が対象 `(E, ψ)` に割り当てる
`Φ.map (Base ψ)⁻¹ (Div ψ)` のことである。 -/
theorem overVal_comp_iso {C' : Type u2} [Category.{v2} C'] {Φ' : MonoidOn.{v, u, w} D}
    (Q : PreFrobenioid C' Φ') {A B E : C'} (ψ : A ⟶ B) (e : B ⟶ E) [IsIso e]
    (hψ : IsIso (Q.Base ψ)) (hbe : IsIso (Q.Base e)) (hc : IsIso (Q.Base (ψ ≫ e))) :
    Φ'.map (@inv _ _ _ _ (Q.Base (ψ ≫ e)) hc) (Q.Div (ψ ≫ e))
      = Φ'.map (@inv _ _ _ _ (Q.Base e) hbe) (Φ'.map (@inv _ _ _ _ (Q.Base ψ) hψ) (Q.Div ψ)) := by
  have hdiv : Q.Div (ψ ≫ e) = Q.Div ψ := by
    rw [Q.Div_comp, show Q.Div e = 0 from isIsometric_of_isIso Q e, map_zero, zero_add,
      show Q.degFr e = 1 from isLinear_of_isIso Q e]
    simp
  have hinv : @inv _ _ _ _ (Q.Base (ψ ≫ e)) hc
      = @inv _ _ _ _ (Q.Base e) hbe ≫ @inv _ _ _ _ (Q.Base ψ) hψ := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [Q.Base_comp, Category.assoc, ← Category.assoc (Q.Base e), IsIso.hom_inv_id,
      Category.id_comp]
    exact @IsIso.hom_inv_id _ _ _ _ (Q.Base ψ) hψ
  rw [hdiv, hinv]
  exact Φ'.map_comp _ _ _

end Abstract

/-! ## ★2. `Λ_k` は底を変えない -/

section Lam

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★**`Λ_k` は底を変えない**。

★`rootBase_mk_spec` を `idxOne` に当て、`rtMap_spec` の底版と突き合わせて
`Base (rtExt B k)`(同型)で割る。

★**型をきれいにして受ける**のが要点 —— `rootBase_mk_spec` をそのまま
`have` で受けると `show` の包みが残り、`rw` が通らない。 -/
theorem rootBase_lamHom (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    (pfRootPre P F).Base (lamHom (F := F) k φ) = P.Base φ := by
  haveI hB : IsIso (P.Base (rtExt P F B k)) := (rtExt_frobType P F B k).2
  have h1 : rootBase (lamHom (F := F) k φ) ≫ P.Base (rtExt P F B k) ≫ P.Base (𝟙 (rtObj P F B k))
      = (P.Base (rtExt P F A k) ≫ P.Base (𝟙 (rtObj P F A k))) ≫ P.Base (rtMap (F := F) k φ) :=
    rootBase_mk_spec (X := (⟨A, k⟩ : PfRootObj P F)) (Y := (⟨B, k⟩ : PfRootObj P F))
      (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)
  rw [P.Base_id, P.Base_id, Category.comp_id, Category.comp_id] at h1
  have h2 := congrArg P.Base (rtMap_spec (F := F) k φ)
  rw [P.Base_comp, P.Base_comp] at h2
  show rootBase (lamHom (F := F) k φ) = P.Base φ
  refine (cancel_mono (P.Base (rtExt P F B k))).mp ?_
  show rootBase (lamHom (F := F) k φ) ≫ P.Base (rtExt P F B k) = P.Base φ ≫ _
  rw [h2]
  exact h1

/-- ★`Λ_k` の底の逆射はもとの底の逆射。 -/
theorem inv_rootBase_lamHom (k : ℕ+) {A B : C} (ψ : A ⟶ B) (hψ : IsIso (P.Base ψ))
    (hb : IsIso ((pfRootPre P F).Base (lamHom (F := F) k ψ))) :
    @inv _ _ _ _ ((pfRootPre P F).Base (lamHom (F := F) k ψ)) hb
      = @inv _ _ _ _ (P.Base ψ) hψ := by
  refine IsIso.inv_eq_of_hom_inv_id ?_
  rw [rootBase_lamHom]
  exact @IsIso.hom_inv_id _ _ _ _ (P.Base ψ) hψ

/-- ★★**`Λ_k` の「後置の値」** —— `k` で割られる。 -/
theorem overVal_lamHom (k : ℕ+) {A B : C} (ψ : A ⟶ B) (hψ : IsIso (P.Base ψ))
    (hb : IsIso ((pfRootPre P F).Base (lamHom (F := F) k ψ))) :
    (Φ.pfOn (phiSharp P)).map
        (@inv _ _ _ _ ((pfRootPre P F).Base (lamHom (F := F) k ψ)) hb)
        ((pfRootPre P F).Div (lamHom (F := F) k ψ))
      = Pf.mk (Φ.map (@inv _ _ _ _ (P.Base ψ) hψ) (P.Div ψ)) k := by
  show Pf.map (Φ.map (@inv _ _ _ _ ((pfRootPre P F).Base (lamHom (F := F) k ψ)) hb))
      ((pfRootPre P F).Div (lamHom (F := F) k ψ)) = _
  rw [rootDiv_lamHom, inv_rootBase_lamHom k ψ hψ hb]
  exact Pf.map_mk _ _ _

end Lam

/-! ## ★3. (iii)(d) 後置の本質的全射性 -/

section EssSurj

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

set_option maxHeartbeats 1000000 in
/-- ★★★**(iii)(d) 後置の本質的全射性** —— `𝒞^pf` 版。

★★**手**(前置の場合と対になる): `X = (A,r)` と `x = a/n` に対し
`X ≅ (A^{(n)}, r·n)` で根を上げ、`𝒞` の `coaPreOverEquiv` の本質的全射性を
`A^{(n)}` で当て、その co-angular pre-step を `Λ_{r·n}` で押し出してから
**根を上げる同型の逆で `X` へ戻す**。
★戻す段は `overVal_comp_iso` が処理する。 -/
theorem pfRoot_coaPreOver_essSurj (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P)
    (X : PfRootObj P F) :
    letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
    (coaPreOverFunctor (pfRootPre P F) X).EssSurj := by
  letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  refine ⟨fun xo => ?_⟩
  obtain ⟨x, rfl⟩ : ∃ x, Opposite.op x = xo := ⟨xo.unop, rfl⟩
  obtain ⟨a, n, rfl⟩ : ∃ (a : Φ.val (P.toElem.obj X.obj).base) (n : ℕ+), x = Pf.mk a n :=
    Pf.inductionOn (p := fun y => ∃ (a : Φ.val (P.toElem.obj X.obj).base) (n : ℕ+),
      y = Pf.mk a n) x (fun m b => ⟨m, b, rfl⟩)
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root n (X.root * n) rfl
  haveI := hXiso
  haveI hbiso : IsIso ((pfRootPre P F).Base eX) :=
    isBaseIsomorphism_of_isIso (pfRootPre P F) eX
  obtain ⟨Z₀, ⟨e₀⟩⟩ := (G.coaPreOverEquiv (rtObj P F X.obj n)).essSurj.mem_essImage
    (Opposite.op (toOrderCat (Φ.map (@inv _ _ _ _ ((pfRootPre P F).Base eX) hbiso)
      (((X.root : ℕ+) : ℕ) • a))))
  haveI hψb : IsIso (P.Base Z₀.hom.hom) := Z₀.hom.property.2.2
  have hdiv0 : Φ.map (@inv _ _ _ _ (P.Base Z₀.hom.hom) hψb) (P.Div Z₀.hom.hom)
      = Φ.map (@inv _ _ _ _ ((pfRootPre P F).Base eX) hbiso) (((X.root : ℕ+) : ℕ) • a) :=
    mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2
      (leOfHom e₀.inv.unop) (leOfHom e₀.hom.unop)
  have hlam : IsPreStep (pfRootPre P F) (lamHom (F := F) (X.root * n) Z₀.hom.hom) :=
    lamHom_isPreStep (X.root * n) Z₀.hom.hom Z₀.hom.property.2
  haveI hlamb : IsIso ((pfRootPre P F).Base (lamHom (F := F) (X.root * n) Z₀.hom.hom)) := hlam.2
  haveI hinvb : IsIso ((pfRootPre P F).Base (inv eX)) :=
    isBaseIsomorphism_of_isIso (pfRootPre P F) (inv eX)
  have hφstep : IsPreStep (pfRootPre P F)
      (lamHom (F := F) (X.root * n) Z₀.hom.hom ≫ inv eX) :=
    IsPreStep.comp (pfRootPre P F) hlam (isPreStep_of_isIso (pfRootPre P F) (inv eX))
  haveI hφb : IsIso ((pfRootPre P F).Base
      (lamHom (F := F) (X.root * n) Z₀.hom.hom ≫ inv eX)) := hφstep.2
  refine ⟨Over.mk (show (⟨(⟨Z₀.left.obj, X.root * n⟩ : PfRootObj P F)⟩ :
      WideSubcategory (coaPreProp (pfRootPre P F))) ⟶ ⟨X⟩ from
    ⟨lamHom (F := F) (X.root * n) Z₀.hom.hom ≫ inv eX,
      pfRoot_isCoAngular hfi _, hφstep⟩), ⟨eqToIso ?_⟩⟩
  refine congrArg Opposite.op (congrArg toOrderCat ?_)
  show (Φ.pfOn (phiSharp P)).map
      (@inv _ _ _ _ ((pfRootPre P F).Base
        (lamHom (F := F) (X.root * n) Z₀.hom.hom ≫ inv eX)) hφb)
      ((pfRootPre P F).Div (lamHom (F := F) (X.root * n) Z₀.hom.hom ≫ inv eX))
    = Pf.mk a n
  rw [overVal_comp_iso (pfRootPre P F) (lamHom (F := F) (X.root * n) Z₀.hom.hom) (inv eX)
      hlamb hinvb hφb,
    overVal_lamHom (X.root * n) Z₀.hom.hom hψb hlamb, hdiv0,
    inv_Base_inv_eq (pfRootPre P F) eX hbiso hinvb]
  show Pf.map (Φ.map ((pfRootPre P F).Base eX))
      (Pf.mk (Φ.map (@inv _ _ _ _ ((pfRootPre P F).Base eX) hbiso)
        (((X.root : ℕ+) : ℕ) • a)) (X.root * n)) = Pf.mk a n
  rw [Pf.map_mk, MonoidOn.map_inv_map Φ _ hbiso]
  exact Pf.mk_smul_cancel a X.root n

end EssSurj

/-! ## ★出典の紐付け(条つき) -/

def pfRoot_coaPreOver_essSurj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58,
    item := "Proposition 3.2, (iii) — (iii)(d) 後置の本質的全射性",
    sectionId := "frdi-prop-3-2" }

end ABC3.Found.FrdI
