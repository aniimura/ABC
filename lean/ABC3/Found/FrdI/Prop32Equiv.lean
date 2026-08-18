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

/-- ★`Φ.map g⁻¹ ∘ Φ.map g = id`。 -/
theorem MonoidOn.map_map_inv (Φ : MonoidOn.{v, u, w} D) {X Y : D} (g : X ⟶ Y) (hg : IsIso g)
    (y : Φ.val Y) : Φ.map (@inv _ _ _ _ g hg) (Φ.map g y) = y := by
  rw [← Φ.map_comp, @IsIso.inv_hom_id _ _ _ _ g hg]
  exact Φ.map_id Y y

/-- ★★**同型に沿った `Φ.map` は `≼` を反映する**。

★充満性で「底に沿って落とした因子の `≼`」を「もとの因子の `≼`」に戻すのに使う。
`Base` が同型なのは **Frobenius 型射**だからである。 -/
theorem mle_of_map_mle (Φ : MonoidOn.{v, u, w} D) {X Y : D} (g : X ⟶ Y) (hg : IsIso g)
    {x y : Φ.val Y} (h : MLe (Φ.map g x) (Φ.map g y)) : MLe x y := by
  obtain ⟨z, hz⟩ := h
  refine ⟨Φ.map (@inv _ _ _ _ g hg) z, ?_⟩
  have h2 := congrArg (Φ.map (@inv _ _ _ _ g hg)) hz
  rw [map_add, MonoidOn.map_map_inv Φ g hg, MonoidOn.map_map_inv Φ g hg] at h2
  exact h2

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

/-! ## ★1b. モノイドの側の descent —— `k` 倍の `≼` から `≼` へ

★★これが充満性の算術の核である。`Pf` の `≼` は分子の `k` 倍でしか降りない
(`Pf.mle_num_of_mle`)ので、そこから `k` を落とす段が要る。 -/

section MonoidDescent

universe w'

/-- ★★**divisorial なモノイドでは `k • a ≼ k • b ⟹ a ≼ b`**。

★**saturated**(`k` 倍が像に入るなら本人も像に入る)と
**integral**(`toGp` が単射)の**両方**が要る。
`Definition 1.1, (i)` の `pre-divisorial` はこの 2 つを含んでいる。 -/
theorem mle_of_nsmul_mle {M : Type w'} [AddCommMonoid M] (hsat : IsSaturatedMonoid M)
    (hint : IsIntegralMonoid M) {a b : M} {k : ℕ} (hk : 0 < k)
    (h : MLe (k • a) (k • b)) : MLe a b := by
  obtain ⟨z, hz⟩ := h
  have hgp : ((k : ℕ) • (toGp M b - toGp M a)) = toGp M z := by
    have h1 := congrArg (toGp M) hz
    rw [toGp_add] at h1
    rw [smul_sub, ← toGp_nsmul, ← toGp_nsmul, ← h1]
    abel
  obtain ⟨w, hw⟩ := hsat (toGp M b - toGp M a) k hk ⟨z, hgp.symm⟩
  refine ⟨w, hint ?_⟩
  rw [toGp_add, hw]
  abel

end MonoidDescent

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

/-! ## ★4. 充満性への下ごしらえ —— 同根の場合の零因子の明示式

★★**同じ 3 脚添字の第 1 脚を共有する 2 本の射は、分母も分子の写像も一致する**
——これが (iii)(d) の充満性の要である(`Prop32Frob.lean` の ★58 の測定)。
下の式で `Θ = Φ.map (Base (rtExt A r)) ∘ Φ.map (Base Z.hom.hom.1)` と
`N = r·r·(deg の第 1 脚)` が**添字だけで決まる**ことが見える。 -/

section SameRoot

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★同根の場合の `rootDiv` の明示式。 -/
theorem rootDiv_mk_sameRoot {A B : C} {r : ℕ+}
    (Z : IdxPf P F (rtObj P F A r) (rtObj P F B r)) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    (pfRootPre P F).Div (show HomRoot P F (⟨A, r⟩ : PfRootObj P F) ⟨B, r⟩ from HomPf.mk Z φ)
      = Pf.mk (Φ.map (P.Base (rtExt P F A r)) (Φ.map (P.Base Z.hom.hom.1) (P.Div φ)))
          (r * r * repRoot Z) := by
  show Pf.divBy (r * r) (Pf.map (Φ.map (P.Base (rtExt P F A r))) (pfDiv (HomPf.mk Z φ))) = _
  rw [pfDiv_mk]
  show Pf.divBy (r * r) (Pf.map (Φ.map (P.Base (rtExt P F A r)))
    (Pf.mk (Φ.map (P.Base Z.hom.hom.1) (P.Div φ)) (repRoot Z))) = _
  rw [Pf.map_mk, Pf.divBy_mk]

/-! ## ★5. (iii)(d) 前置の充満性(同根の場合)

★★**手**(4 段):

1. `rtRootIso` で根 `r` を `r·r` へ上げ、`exists_rep3_span_isotropic` で
   **2 射を同じ 3 脚添字の上に**、しかも**第 1 脚が isotropic** になるように揃える。
   ★isotropic にしておくと `𝒞` の側で **pre-step が自動的に co-angular** になる
   (`isCoAngular_of_isotropic_dom`)——`𝒞` の (iii)(d) に渡すのにこれが要る。
2. 因子は**分母も分子の写像も一致する**(`rootDiv_rtRootIso_mk`)ので、
   `Pf` の `≼` から代表元の `≼` へ降りる(`mle_div_of_mle_pf`)。
3. `𝒞` の (iii)(d) の**充満性**を第 1 脚の対象で当てて `m₀` を得る。
4. `compRoot_mk3` で `𝒞^pf` へ戻す。 -/

-- ★`exists_idx3_isotropic`(3 脚添字を isotropic な代表へ押し上げる)は
-- `Prop32Frob.lean` に既にある。ここではそれを使う。

/-- ★★**始域を共有する 2 射を、isotropic な代表を持つ同じ 3 脚添字へ**。 -/
theorem exists_rep3_span_isotropic (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (f : HomPf P F A B) (g : HomPf P F A E) :
    ∃ (V : IdxPf3 P F A B E) (φ : V.right.obj.1 ⟶ V.right.obj.2.1)
      (ψ : V.right.obj.1 ⟶ V.right.obj.2.2),
      IsIsotropic P V.right.obj.1 ∧
      f = HomPf.mk ((idx12 P F A B E).obj V) φ ∧
      g = HomPf.mk ((idx13 P F A B E).obj V) ψ := by
  obtain ⟨V₀, φ₀, ψ₀, hf₀, hg₀⟩ := exists_rep3_span (P := P) (F := F) f g
  obtain ⟨V, u, hV⟩ := exists_idx3_isotropic (F := F) hfi V₀
  refine ⟨V, idxTransport P F ((idx12 P F A B E).map u) φ₀,
    idxTransport P F ((idx13 P F A B E).map u) ψ₀, hV, ?_, ?_⟩
  · rw [hf₀]
    exact (HomPf.mk_map ((idx12 P F A B E).map u) φ₀).symm
  · rw [hg₀]
    exact (HomPf.mk_map ((idx13 P F A B E).map u) ψ₀).symm

set_option maxHeartbeats 1000000 in
/-- ★★**充満性の算術** —— `Pf` の `≼` を代表元の `≼` に降ろす。

★3 段: `Pf.mle_num_of_mle`(分子の `k` 倍まで)→ `mle_of_map_mle` を 2 回
(2 つの底はいずれも Frobenius 型射の底なので同型)→ `mle_of_nsmul_mle`
(divisorial の saturated ＋ integral)。 -/
theorem mle_div_of_mle_pf {A B₁ B₂ : C} {r : ℕ+}
    (V : IdxPf3 P F (rtObj P F A (r*r)) (rtObj P F B₁ (r*r)) (rtObj P F B₂ (r*r)))
    (f : V.right.obj.1 ⟶ V.right.obj.2.1) (g : V.right.obj.1 ⟶ V.right.obj.2.2)
    (hle : MLe
      (Pf.mk (Φ.map (P.Base (rtExt P F A (r * r))) (Φ.map (P.Base V.hom.hom.1) (P.Div f)))
        (r * r * (P.degFr V.hom.hom.1 * r)))
      (Pf.mk (Φ.map (P.Base (rtExt P F A (r * r))) (Φ.map (P.Base V.hom.hom.1) (P.Div g)))
        (r * r * (P.degFr V.hom.hom.1 * r)))) :
    MLe (P.Div f) (P.Div g) := by
  haveI h1 : IsIso (P.Base (rtExt P F A (r*r))) := (rtExt_frobType P F A (r*r)).2
  haveI h2 : IsIso (P.Base V.hom.hom.1) := V.hom.property.1.2
  obtain ⟨k, hk⟩ := Pf.mle_num_of_mle hle
  rw [← map_nsmul, ← map_nsmul] at hk
  have hk2 := mle_of_map_mle Φ (P.Base (rtExt P F A (r*r))) h1 hk
  have hk2' : MLe (Φ.map (P.Base V.hom.hom.1) (((k : ℕ+) : ℕ) • P.Div f))
      (Φ.map (P.Base V.hom.hom.1) (((k : ℕ+) : ℕ) • P.Div g)) := by
    rw [map_nsmul, map_nsmul]
    exact hk2
  exact mle_of_nsmul_mle (P.divisorial _).1.2.1 (P.divisorial _).1.1 k.2
    (mle_of_map_mle Φ (P.Base V.hom.hom.1) h2 hk2')

set_option maxHeartbeats 2000000 in
/-- ★★★**(iii)(d) 前置の充満性(同根の場合)**。 -/
theorem pfRoot_coaPreUnder_full_sameRoot (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P)
    {A B₁ B₂ : C} {r : ℕ+}
    (φ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B₁, r⟩) (ψ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B₂, r⟩)
    (hφs : IsPreStep (pfRootPre P F) φ) (hψs : IsPreStep (pfRootPre P F) ψ)
    (hle : MLe ((pfRootPre P F).Div φ) ((pfRootPre P F).Div ψ)) :
    ∃ m : (⟨B₁, r⟩ : PfRootObj P F) ⟶ ⟨B₂, r⟩,
      IsPreStep (pfRootPre P F) m ∧ φ ≫ m = ψ := by
  obtain ⟨V, f, g, hViso, hf, hg⟩ := exists_rep3_span_isotropic (F := F) hfi
    ((rtRootIso P F A B₁ (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ)
    ((rtRootIso P F A B₂ (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv ψ)
  have hφmk : φ = (rtRootIso P F A B₁ (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx12 P F _ _ _).obj V) f) := by
    rw [← hf, Iso.inv_hom_id_apply]
  have hψmk : ψ = (rtRootIso P F A B₂ (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx13 P F _ _ _).obj V) g) := by
    rw [← hg, Iso.inv_hom_id_apply]
  have hfs : IsPreStep P f := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B₁, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F B₁ (show r * r = r * r from rfl))
        (rtLift_frobType P F B₁ _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx12 P F _ _ _).obj V)) f).mp ?_
    rw [hφmk, rtRootIso_hom_mk] at hφs
    exact hφs
  have hgs : IsPreStep P g := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B₂, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F B₂ (show r * r = r * r from rfl))
        (rtLift_frobType P F B₂ _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx13 P F _ _ _).obj V)) g).mp ?_
    rw [hψmk, rtRootIso_hom_mk] at hψs
    exact hψs
  have hfc : IsCoAngular P f := isCoAngular_of_isotropic_dom (P := P) F hViso f
  have hgc : IsCoAngular P g := isCoAngular_of_isotropic_dom (P := P) F hViso g
  have hdivφ : (pfRootPre P F).Div φ
      = Pf.mk (Φ.map (P.Base (rtExt P F A (r * r))) (Φ.map (P.Base V.hom.hom.1) (P.Div f)))
          (r * r * (P.degFr V.hom.hom.1 * r)) := by
    rw [hφmk]
    exact rootDiv_rtRootIso_mk (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B₁, r⟩ : PfRootObj P F))
      (show r * r = r * r from rfl) (show r * r = r * r from rfl) ((idx12 P F _ _ _).obj V) f
  have hdivψ : (pfRootPre P F).Div ψ
      = Pf.mk (Φ.map (P.Base (rtExt P F A (r * r))) (Φ.map (P.Base V.hom.hom.1) (P.Div g)))
          (r * r * (P.degFr V.hom.hom.1 * r)) := by
    rw [hψmk]
    exact rootDiv_rtRootIso_mk (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B₂, r⟩ : PfRootObj P F))
      (show r * r = r * r from rfl) (show r * r = r * r from rfl) ((idx13 P F _ _ _).obj V) g
  rw [hdivφ, hdivψ] at hle
  have hmle : MLe (P.Div f) (P.Div g) := mle_div_of_mle_pf V f g hle
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv V.right.obj.1
  obtain ⟨uu, -⟩ := (coaPreUnderFunctor P V.right.obj.1).map_surjective
    (show (coaPreUnderFunctor P V.right.obj.1).obj
        (Under.mk (show (⟨V.right.obj.1⟩ : WideSubcategory (coaPreProp P))
          ⟶ ⟨V.right.obj.2.1⟩ from ⟨f, hfc, hfs⟩))
      ⟶ (coaPreUnderFunctor P V.right.obj.1).obj
        (Under.mk (show (⟨V.right.obj.1⟩ : WideSubcategory (coaPreProp P))
          ⟶ ⟨V.right.obj.2.2⟩ from ⟨g, hgc, hgs⟩))
      from homOfLE hmle)
  have hcomp : f ≫ uu.right.hom = g := congrArg InducedWideCategory.Hom.hom (Under.w uu)
  refine ⟨(rtRootIso P F B₁ B₂ (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx23 P F _ _ _).obj V) uu.right.hom), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    exact (isPreStep_mk_iff (X := (⟨B₁, r⟩ : PfRootObj P F)) (Y := (⟨B₂, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F B₁ (show r * r = r * r from rfl))
        (rtLift_frobType P F B₁ _) (rtLift P F B₂ (show r * r = r * r from rfl))
        (rtLift_frobType P F B₂ _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx23 P F _ _ _).obj V)) uu.right.hom).mpr uu.right.property.2
  · rw [hφmk, hψmk]
    refine Eq.trans ?_ (congrArg (fun t : V.right.obj.1 ⟶ V.right.obj.2.2 =>
      (rtRootIso P F A B₂ (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).hom (HomPf.mk ((idx13 P F _ _ _).obj V) t)) hcomp)
    exact compRoot_mk3 (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B₁, r⟩ : PfRootObj P F))
      (Z := (⟨B₂, r⟩ : PfRootObj P F)) V.hom.hom.1 V.hom.hom.2.1 V.hom.hom.2.2
      V.hom.property.1 V.hom.property.2.1 V.hom.property.2.2.1
      V.hom.property.2.2.2.1 V.hom.property.2.2.2.2 f uu.right.hom

/-! ## ★6. (iii)(d) 前置の充満性(一般) -/

-- ★`MLe.map`(`≼` は加法準同型で保たれる)は `Prop110.lean` に既にある。

/-- ★同型で挟むと `Div` は底に沿って写るだけ。 -/
theorem rootDiv_conj {X X' Y Y' : PfRootObj P F} (e : X' ⟶ X) (h : X ⟶ Y) (e' : Y ⟶ Y')
    [IsIso e] [IsIso e'] :
    (pfRootPre P F).Div (e ≫ h ≫ e')
      = (Φ.pfOn (phiSharp P)).map ((pfRootPre P F).Base e) ((pfRootPre P F).Div h) := by
  rw [(pfRootPre P F).Div_comp, (pfRootPre P F).Div_comp,
    show (pfRootPre P F).Div e' = 0 from isIsometric_of_isIso (pfRootPre P F) e',
    show (pfRootPre P F).Div e = 0 from isIsometric_of_isIso (pfRootPre P F) e,
    show (pfRootPre P F).degFr e' = 1 from isLinear_of_isIso (pfRootPre P F) e']
  simp

set_option maxHeartbeats 1000000 in
/-- ★★★**(iii)(d) 前置の充満性(射の形)** —— 根を揃えて同根版に落とす。

★`pfRoot_faithfulUpToUnits` と同じ共役の定型:
`X`・`Y₁`・`Y₂` の根をすべて `X.root · Y₁.root · Y₂.root` へ上げ、
`inv eX ≫ − ≫ e_i` で移してから同根版を当て、`e₁ ≫ − ≫ inv e₂` で戻す。 -/
theorem pfRoot_coaPreUnder_full_hom (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P)
    {X Y₁ Y₂ : PfRootObj P F} (φ : X ⟶ Y₁) (ψ : X ⟶ Y₂)
    (hφs : IsPreStep (pfRootPre P F) φ) (hψs : IsPreStep (pfRootPre P F) ψ)
    (hle : MLe ((pfRootPre P F).Div φ) ((pfRootPre P F).Div ψ)) :
    ∃ m : Y₁ ⟶ Y₂, IsPreStep (pfRootPre P F) m ∧ φ ≫ m = ψ := by
  obtain ⟨eX, hX⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root (Y₁.root * Y₂.root)
    (X.root * Y₁.root * Y₂.root) (by ac_rfl)
  obtain ⟨e1, h1⟩ := pfRoot_exists_iso_root (F := F) Y₁.obj Y₁.root (X.root * Y₂.root)
    (X.root * Y₁.root * Y₂.root) (by ac_rfl)
  obtain ⟨e2, h2⟩ := pfRoot_exists_iso_root (F := F) Y₂.obj Y₂.root (X.root * Y₁.root)
    (X.root * Y₁.root * Y₂.root) (by ac_rfl)
  haveI := hX; haveI := h1; haveI := h2
  have hφ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ e1) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) e1))
  have hψ' : IsPreStep (pfRootPre P F) (inv eX ≫ ψ ≫ e2) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hψs (isPreStep_of_isIso (pfRootPre P F) e2))
  have hle' : MLe ((pfRootPre P F).Div (inv eX ≫ φ ≫ e1))
      ((pfRootPre P F).Div (inv eX ≫ ψ ≫ e2)) := by
    rw [rootDiv_conj, rootDiv_conj]
    exact MLe.map _ hle
  obtain ⟨m', hm's, hm'e⟩ := pfRoot_coaPreUnder_full_sameRoot hfi G
    (inv eX ≫ φ ≫ e1) (inv eX ≫ ψ ≫ e2) hφ' hψ' hle'
  refine ⟨e1 ≫ m' ≫ inv e2, ?_, ?_⟩
  · exact IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) e1)
      (IsPreStep.comp (pfRootPre P F) hm's (isPreStep_of_isIso (pfRootPre P F) (inv e2)))
  · have h3 : eX ≫ (inv eX ≫ φ ≫ e1) ≫ m' = eX ≫ (inv eX ≫ ψ ≫ e2) := by rw [hm'e]
    simp only [Category.assoc, IsIso.hom_inv_id_assoc] at h3
    have h4 : (φ ≫ e1 ≫ m') ≫ inv e2 = (ψ ≫ e2) ≫ inv e2 :=
      congrArg (fun t => t ≫ inv e2) h3
    simpa only [Category.assoc, IsIso.hom_inv_id, Category.comp_id] using h4

set_option maxHeartbeats 1000000 in
/-- ★★★★**(iii)(d) 前置の充満性**。 -/
theorem pfRoot_coaPreUnder_full (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P)
    (X : PfRootObj P F) :
    letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
    (coaPreUnderFunctor (pfRootPre P F) X).Full := by
  letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
  refine ⟨fun {Z W} h => ?_⟩
  obtain ⟨m, hms, hme⟩ := pfRoot_coaPreUnder_full_hom hfi G Z.hom.hom W.hom.hom
    Z.hom.property.2 W.hom.property.2 (leOfHom h)
  refine ⟨Under.homMk (show (⟨Z.right.obj⟩ : WideSubcategory (coaPreProp (pfRootPre P F)))
      ⟶ ⟨W.right.obj⟩ from ⟨m, pfRoot_isCoAngular hfi m, hms⟩)
    (WideSubcategory.hom_ext _ hme), Subsingleton.elim _ _⟩

end SameRoot

/-! ## ★出典の紐付け(条つき) -/

def pfRoot_coaPreOver_essSurj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58,
    item := "Proposition 3.2, (iii) — (iii)(d) 後置の本質的全射性",
    sectionId := "frdi-prop-3-2" }

end ABC3.Found.FrdI
