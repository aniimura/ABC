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
**このファイルは残る 3 条**(後置の本質的全射性・前置の充満性・後置の充満性)を埋め、
★★★**`pfRoot_frobenioid`(`𝒞^pf` は Frobenioid)を組み立てる**。

## ★★★2026-08-18 —— (iii)(d) の 6 条が揃った

| 条 | 宣言 | 手 |
|---|---|---|
| 前置・忠実 | `coaPreUnder_faithful` | pre-Frobenioid 一般(全射性) |
| 後置・忠実 | `coaPreOver_faithful` | pre-Frobenioid 一般(pre-step が mono) |
| 前置・本質的全射 | `pfRoot_coaPreUnder_essSurj` | 根を上げて `Λ_k` で押し出す |
| 後置・本質的全射 | `pfRoot_coaPreOver_essSurj` | 同上 ＋ `overVal_comp_iso` で戻す |
| 前置・充満 | `pfRoot_coaPreUnder_full` | ★3 脚添字(スパン)＋ `mle_div_of_mle_pf` |
| 後置・充満 | ★`pfRoot_coaPreOver_full` | ★3 脚添字(**余スパン**)＋ `mle_of_mle_pf_map` |

★★充満性 2 本の共通の急所は「**同じ 3 脚添字の上では、2 本の射の因子が
同じ写像・同じ分母で書ける**」ことである(`rootDiv_rtRootIso_mk` /
★`overVal_rtRootIso_mk`)。そこまで来れば `Pf` の `≼` は代表元の `≼` に降りる。

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

/-- ★★**後置の値は始域側の同型を無視する** —— `overVal (e ≫ ψ) = overVal ψ`。

★「後置」はスライス `(𝒞^coa-pre)_A` の値なので、**終域側の同型でだけ動く**
(`overVal_comp_iso`)。始域側は素通りする。
★根を揃える共役 `inv eX ≫ ζ ≫ eA` の `≼` を移すのにこの 2 本が要る。 -/
theorem overVal_iso_comp {C' : Type u2} [Category.{v2} C'] {Φ' : MonoidOn.{v, u, w} D}
    (Q : PreFrobenioid C' Φ') {X X' Y : C'} (e : X' ⟶ X) (ψ : X ⟶ Y) [IsIso e]
    (hψ : IsIso (Q.Base ψ)) (hbe : IsIso (Q.Base e)) (hc : IsIso (Q.Base (e ≫ ψ))) :
    Φ'.map (@inv _ _ _ _ (Q.Base (e ≫ ψ)) hc) (Q.Div (e ≫ ψ))
      = Φ'.map (@inv _ _ _ _ (Q.Base ψ) hψ) (Q.Div ψ) := by
  have hdiv : Q.Div (e ≫ ψ) = Φ'.map (Q.Base e) (Q.Div ψ) := by
    rw [Q.Div_comp, show Q.Div e = 0 from isIsometric_of_isIso Q e]
    simp
  have hinv : @inv _ _ _ _ (Q.Base (e ≫ ψ)) hc
      = @inv _ _ _ _ (Q.Base ψ) hψ ≫ @inv _ _ _ _ (Q.Base e) hbe := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [Q.Base_comp, Category.assoc, ← Category.assoc (Q.Base ψ), IsIso.hom_inv_id,
      Category.id_comp]
    exact @IsIso.hom_inv_id _ _ _ _ (Q.Base e) hbe
  rw [hdiv, hinv, Φ'.map_comp]
  exact congrArg (Φ'.map (@inv _ _ _ _ (Q.Base ψ) hψ))
    (MonoidOn.map_map_inv Φ' (Q.Base e) hbe (Q.Div ψ))

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

/-! ## ★7. (iii)(d) 後置の充満性への下ごしらえ —— 「後置の値」の明示式

★★前置(`rootDiv_rtRootIso_mk`)の**後置版**。要点は

  `overVal c = Pf.mk (Φ.map u (Φ.map (Base φ)⁻¹ (Div φ))) N`,
  `u = Base(rtExt Y.obj tB) ≫ Base W.2`,  `N = X.root · Y.root · (degFr W.1 · e)`

で、★★**`u` は第 2 脚(＝余スパンの共有側)だけで決まる**こと。
余スパンの 2 本には同じ 3 脚添字から同じ `u` が出るので、前置と同じ降下
(`Pf.mle_num_of_mle` → `mle_of_map_mle` → `mle_of_nsmul_mle`)が回る。
★分母は `degFr W.1` が第 1 脚に依るが、3 脚添字では**3 本の脚の次数が等しい**
(`V.hom.property.2.2.2.1`)ので、これも一致する。 -/

/-- ★型の罠を避けるための小さな抽象補題 —— `(g₁ ≫ g₂ ≫ f) ≫ (f⁻¹ ≫ (g₁≫g₂)⁻¹) = 𝟙`。

★`reassoc_of%` は `WideSubcategory` 由来の項では
「target expression is not type-correct」で失敗する。抽象化して逃げる。 -/
theorem comp3_inv_inv_id {D' : Type*} [Category D'] {U V T S : D'}
    (g₁ : U ⟶ V) (g₂ : V ⟶ T) (f : T ⟶ S) [IsIso g₁] [IsIso g₂] [IsIso f] :
    (g₁ ≫ g₂ ≫ f) ≫ (inv f ≫ inv (g₁ ≫ g₂)) = 𝟙 U := by
  simp

/-- ★`rootDiv_rtRootIso_mk` の `Φ.map` を 1 本にまとめた形。 -/
theorem rootDiv_rtRootIso_mk' {X Y : PfRootObj P F} {e tA tB : ℕ+}
    (hA : tA = e * Y.root) (hB : tB = e * X.root)
    (W : IdxPf P F (rtObj P F X.obj tA) (rtObj P F Y.obj tB))
    (φ : W.right.obj.1 ⟶ W.right.obj.2) :
    (pfRootPre P F).Div (show HomRoot P F X Y from
        (rtRootIso P F X.obj Y.obj hA hB).hom (HomPf.mk W φ))
      = Pf.mk (Φ.map (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1) (P.Div φ))
          (X.root * Y.root * (P.degFr W.hom.hom.1 * e)) := by
  refine Eq.trans (rootDiv_rtRootIso_mk hA hB W φ) ?_
  exact congrArg (fun t => Pf.mk t (X.root * Y.root * (P.degFr W.hom.hom.1 * e)))
    (Φ.map_comp _ _ _).symm

set_option maxHeartbeats 1000000 in
/-- ★★★**根の取り替えを通した「後置の値」の明示式**。

★★射を素の変数 `cc` で受けるのが要点(`WideSubcategory` の型の罠を避ける)。 -/
theorem overVal_rtRootIso_mk {X Y : PfRootObj P F} {e tA tB : ℕ+}
    (hA : tA = e * Y.root) (hB : tB = e * X.root)
    (W : IdxPf P F (rtObj P F X.obj tA) (rtObj P F Y.obj tB))
    (φ : W.right.obj.1 ⟶ W.right.obj.2)
    (cc : X ⟶ Y) (hcc : cc = (rtRootIso P F X.obj Y.obj hA hB).hom (HomPf.mk W φ))
    (hφb : IsIso (P.Base φ)) (hcb : IsIso ((pfRootPre P F).Base cc)) :
    (Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base cc) hcb)
        ((pfRootPre P F).Div cc)
      = Pf.mk (Φ.map (P.Base (rtExt P F Y.obj tB) ≫ P.Base W.hom.hom.2)
          (Φ.map (@inv _ _ _ _ (P.Base φ) hφb) (P.Div φ)))
          (X.root * Y.root * (P.degFr W.hom.hom.1 * e)) := by
  haveI hw1 : IsIso (P.Base (rtExt P F X.obj tA)) := (rtExt_frobType P F X.obj tA).2
  haveI hw2 : IsIso (P.Base W.hom.hom.1) := W.hom.property.1.2
  haveI hwc : IsIso (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1) :=
    IsIso.comp_isIso' hw1 hw2
  haveI := hφb
  have hcu' : (pfRootPre P F).Base cc ≫ P.Base (rtExt P F Y.obj tB) ≫ P.Base W.hom.hom.2
      = (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1) ≫ P.Base φ := by
    rw [hcc]
    exact rootBase_rtRootIso_mk_spec hA hB W φ
  have hdiv : (pfRootPre P F).Div cc
      = Pf.mk (Φ.map (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1) (P.Div φ))
          (X.root * Y.root * (P.degFr W.hom.hom.1 * e)) := by
    rw [hcc]
    exact rootDiv_rtRootIso_mk' hA hB W φ
  have hinv : @inv _ _ _ _ ((pfRootPre P F).Base cc) hcb
      = (P.Base (rtExt P F Y.obj tB) ≫ P.Base W.hom.hom.2)
        ≫ inv (P.Base φ) ≫ inv (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1) := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    have key := congrArg
      (fun t => t ≫ (inv (P.Base φ) ≫ inv (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1)))
      hcu'
    simp only [Category.assoc] at key ⊢
    refine key.trans ?_
    exact @comp3_inv_inv_id _ _ _ _ _ _ (P.Base (rtExt P F X.obj tA)) (P.Base W.hom.hom.1)
      (P.Base φ) hw1 hw2 hφb
  rw [hinv, hdiv]
  show Pf.map (Φ.map _) (Pf.mk _ _) = _
  rw [Pf.map_mk]
  congr 1
  refine Eq.trans (Φ.map_comp (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1)
    ((P.Base (rtExt P F Y.obj tB) ≫ P.Base W.hom.hom.2) ≫ inv (P.Base φ)
      ≫ inv (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1)) (P.Div φ)).symm ?_
  refine Eq.trans ?_ (Φ.map_comp (@inv _ _ _ _ (P.Base φ) hφb)
    (P.Base (rtExt P F Y.obj tB) ≫ P.Base W.hom.hom.2) (P.Div φ))
  congr 1
  simp

/-! ## ★8. (iii)(d) 後置の充満性

★★**手**(前置と同じ 4 段。違うのは「余スパン」であること):

1. `rtRootIso` で根 `r` を `r·r` へ上げ、`exists_rep3_cospan_isotropic` で
   **2 射を同じ 3 脚添字の上に**、しかも**第 1・第 2 脚が isotropic** になるように揃える。
   ★余スパンでは `𝒞` に渡す 2 射の**始域が違う**ので、
   isotropic にすべき脚が 2 本になる(`exists_idx3_isotropic12`)。
2. 「後置の値」は**同じ `u` と同じ分母**を持つ(`overVal_rtRootIso_mk`)ので、
   `Pf` の `≼` から代表元の `≼` へ降りる(`mle_of_mle_pf_map`)。
   ★分母が一致するのは 3 脚添字の**3 本の脚の次数が等しい**からである。
3. `𝒞` の (iii)(d) の**後置の充満性**を第 3 脚の対象で当てて `m₀` を得る。
4. `compRoot_mk3` で `𝒞^pf` へ戻す。 -/

/-- ★★3 脚添字の**第 1・第 2 脚の先を同時に isotropic** にできる。 -/
theorem exists_idx3_isotropic12 (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (V : IdxPf3 P F A B E) :
    ∃ (W : IdxPf3 P F A B E) (_ : V ⟶ W),
      IsIsotropic P W.right.obj.1 ∧ IsIsotropic P W.right.obj.2.1 := by
  obtain ⟨W₁, u₁, h₁⟩ := exists_idx3_isotropic (F := F) hfi V
  obtain ⟨W₂, u₂, h₂⟩ := exists_idx3_isotropic2 (F := F) hfi W₁
  exact ⟨W₂, u₁ ≫ u₂, F.isotropicClosed u₂.right.hom.1 h₁, h₂⟩

/-- ★★**終域を共有する 2 射を、第 1・第 2 脚が isotropic な同じ 3 脚添字へ**。 -/
theorem exists_rep3_cospan_isotropic (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (f : HomPf P F A E) (g : HomPf P F B E) :
    ∃ (V : IdxPf3 P F A B E) (θ : V.right.obj.1 ⟶ V.right.obj.2.2)
      (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2),
      IsIsotropic P V.right.obj.1 ∧ IsIsotropic P V.right.obj.2.1 ∧
      f = HomPf.mk ((idx13 P F A B E).obj V) θ ∧
      g = HomPf.mk ((idx23 P F A B E).obj V) ψ := by
  obtain ⟨V₀, θ₀, ψ₀, hf₀, hg₀⟩ := exists_rep3_cospan (P := P) (F := F) f g
  obtain ⟨V, u, hV1, hV2⟩ := exists_idx3_isotropic12 (F := F) hfi V₀
  refine ⟨V, idxTransport P F ((idx13 P F A B E).map u) θ₀,
    idxTransport P F ((idx23 P F A B E).map u) ψ₀, hV1, hV2, ?_, ?_⟩
  · rw [hf₀]
    exact (HomPf.mk_map ((idx13 P F A B E).map u) θ₀).symm
  · rw [hg₀]
    exact (HomPf.mk_map ((idx23 P F A B E).map u) ψ₀).symm

/-- ★★**後置の充満性の算術** —— `Pf` の `≼` を代表元の `≼` に降ろす。

★前置版(`mle_div_of_mle_pf`)と同じ 3 段だが、
**分子の写像が 1 本の同型 `u` にまとまっている**ぶん短い。
★`Φ` は反変なので `Φ.map u : Φ.val T → Φ.val S` の向きに注意。 -/
theorem mle_of_mle_pf_map (Q : PreFrobenioid C Φ) {S T : D} (u : S ⟶ T) (hu : IsIso u)
    (N : ℕ+) {a b : Φ.val T}
    (hle : MLe (Pf.mk (Φ.map u a) N) (Pf.mk (Φ.map u b) N)) : MLe a b := by
  obtain ⟨k, hk⟩ := Pf.mle_num_of_mle hle
  rw [← map_nsmul, ← map_nsmul] at hk
  exact mle_of_nsmul_mle (Q.divisorial T).1.2.1 (Q.divisorial T).1.1 k.2
    (mle_of_map_mle Φ u hu hk)

set_option maxHeartbeats 2000000 in
/-- ★★★**(iii)(d) 後置の充満性(同根の場合)**。 -/
theorem pfRoot_coaPreOver_full_sameRoot (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P)
    {A B E : C} {r : ℕ+}
    (ζ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) (ω : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩)
    (hζs : IsPreStep (pfRootPre P F) ζ) (hωs : IsPreStep (pfRootPre P F) ω)
    (hle : MLe ((Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base ω) hωs.2)
        ((pfRootPre P F).Div ω))
      ((Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base ζ) hζs.2)
        ((pfRootPre P F).Div ζ))) :
    ∃ m : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩,
      IsPreStep (pfRootPre P F) m ∧ m ≫ ω = ζ := by
  obtain ⟨V, θ, ψ, hV1, hV2, hf, hg⟩ := exists_rep3_cospan_isotropic (F := F) hfi
    ((rtRootIso P F A E (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv ζ)
    ((rtRootIso P F B E (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv ω)
  have hζmk : ζ = (rtRootIso P F A E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx13 P F _ _ _).obj V) θ) := by
    rw [← hf, Iso.inv_hom_id_apply]
  have hωmk : ω = (rtRootIso P F B E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx23 P F _ _ _).obj V) ψ) := by
    rw [← hg, Iso.inv_hom_id_apply]
  have hθs : IsPreStep P θ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨E, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx13 P F _ _ _).obj V)) θ).mp ?_
    rw [hζmk, rtRootIso_hom_mk] at hζs
    exact hζs
  have hψs : IsPreStep P ψ := by
    refine (isPreStep_mk_iff (X := (⟨B, r⟩ : PfRootObj P F)) (Y := (⟨E, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F B (show r * r = r * r from rfl))
        (rtLift_frobType P F B _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx23 P F _ _ _).obj V)) ψ).mp ?_
    rw [hωmk, rtRootIso_hom_mk] at hωs
    exact hωs
  have hθc : IsCoAngular P θ := isCoAngular_of_isotropic_dom (P := P) F hV1 θ
  have hψc : IsCoAngular P ψ := isCoAngular_of_isotropic_dom (P := P) F hV2 ψ
  -- ★2 本の「後置の値」——**同じ `u`・同じ分母**で書ける
  have hoζ' : (Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base ζ) hζs.2)
        ((pfRootPre P F).Div ζ)
      = Pf.mk (Φ.map (P.Base (rtExt P F E (r * r)) ≫ P.Base V.hom.hom.2.2)
          (Φ.map (@inv _ _ _ _ (P.Base θ) hθs.2) (P.Div θ)))
          (r * r * (P.degFr V.hom.hom.1 * r)) :=
    overVal_rtRootIso_mk (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨E, r⟩ : PfRootObj P F)) (show r * r = r * r from rfl) (show r * r = r * r from rfl)
      ((idx13 P F _ _ _).obj V) θ ζ hζmk hθs.2 hζs.2
  have hoω' : (Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base ω) hωs.2)
        ((pfRootPre P F).Div ω)
      = Pf.mk (Φ.map (P.Base (rtExt P F E (r * r)) ≫ P.Base V.hom.hom.2.2)
          (Φ.map (@inv _ _ _ _ (P.Base ψ) hψs.2) (P.Div ψ)))
          (r * r * (P.degFr V.hom.hom.1 * r)) := by
    rw [V.hom.property.2.2.2.1]
    exact overVal_rtRootIso_mk (X := (⟨B, r⟩ : PfRootObj P F))
      (Y := (⟨E, r⟩ : PfRootObj P F)) (show r * r = r * r from rfl) (show r * r = r * r from rfl)
      ((idx23 P F _ _ _).obj V) ψ ω hωmk hψs.2 hωs.2
  rw [hoζ', hoω'] at hle
  haveI hu : IsIso (P.Base (rtExt P F E (r * r)) ≫ P.Base V.hom.hom.2.2) :=
    IsIso.comp_isIso' (rtExt_frobType P F E (r * r)).2 V.hom.property.2.2.1.2
  have hmle := mle_of_mle_pf_map P _ hu _ hle
  -- ★`𝒞` の (iii)(d) 後置の充満性を第 3 脚の対象で当てる
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv V.right.obj.2.2
  obtain ⟨uu, -⟩ := (coaPreOverFunctor P V.right.obj.2.2).map_surjective
    (show (coaPreOverFunctor P V.right.obj.2.2).obj
        (Over.mk (show (⟨V.right.obj.1⟩ : WideSubcategory (coaPreProp P))
          ⟶ ⟨V.right.obj.2.2⟩ from ⟨θ, hθc, hθs⟩))
      ⟶ (coaPreOverFunctor P V.right.obj.2.2).obj
        (Over.mk (show (⟨V.right.obj.2.1⟩ : WideSubcategory (coaPreProp P))
          ⟶ ⟨V.right.obj.2.2⟩ from ⟨ψ, hψc, hψs⟩))
      from (homOfLE hmle).op)
  have hcomp : uu.left.hom ≫ ψ = θ := congrArg InducedWideCategory.Hom.hom (Over.w uu)
  refine ⟨(rtRootIso P F A B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx12 P F _ _ _).obj V) uu.left.hom), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    exact (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F B (show r * r = r * r from rfl))
        (rtLift_frobType P F B _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx12 P F _ _ _).obj V)) uu.left.hom).mpr uu.left.property.2
  · rw [hζmk, hωmk]
    refine Eq.trans ?_ (congrArg (fun t : V.right.obj.1 ⟶ V.right.obj.2.2 =>
      (rtRootIso P F A E (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).hom (HomPf.mk ((idx13 P F _ _ _).obj V) t)) hcomp)
    exact compRoot_mk3 (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      (Z := (⟨E, r⟩ : PfRootObj P F)) V.hom.hom.1 V.hom.hom.2.1 V.hom.hom.2.2
      V.hom.property.1 V.hom.property.2.1 V.hom.property.2.2.1
      V.hom.property.2.2.2.1 V.hom.property.2.2.2.2 uu.left.hom ψ

set_option maxHeartbeats 2000000 in
/-- ★★★**(iii)(d) 後置の充満性(射の形)** —— 根を揃えて同根版に落とす。

★前置版(`pfRoot_coaPreUnder_full_hom`)と対になる共役の定型。
★「後置の値」は**始域側の同型を無視し、終域側の同型でだけ動く**
(`overVal_iso_comp` / `overVal_comp_iso`)ので、`≼` は `MLe.map` で移せる。 -/
theorem pfRoot_coaPreOver_full_hom (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P)
    {X Y A : PfRootObj P F} (ζ : X ⟶ A) (ω : Y ⟶ A)
    (hζs : IsPreStep (pfRootPre P F) ζ) (hωs : IsPreStep (pfRootPre P F) ω)
    (hle : MLe ((Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base ω) hωs.2)
        ((pfRootPre P F).Div ω))
      ((Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base ζ) hζs.2)
        ((pfRootPre P F).Div ζ))) :
    ∃ m : X ⟶ Y, IsPreStep (pfRootPre P F) m ∧ m ≫ ω = ζ := by
  obtain ⟨eX, hX⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root (Y.root * A.root)
    (X.root * Y.root * A.root) (by ac_rfl)
  obtain ⟨eY, hY⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root (X.root * A.root)
    (X.root * Y.root * A.root) (by ac_rfl)
  obtain ⟨eA, hA⟩ := pfRoot_exists_iso_root (F := F) A.obj A.root (X.root * Y.root)
    (X.root * Y.root * A.root) (by ac_rfl)
  haveI := hX; haveI := hY; haveI := hA
  have hζ' : IsPreStep (pfRootPre P F) (inv eX ≫ ζ ≫ eA) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hζs (isPreStep_of_isIso (pfRootPre P F) eA))
  have hω' : IsPreStep (pfRootPre P F) (inv eY ≫ ω ≫ eA) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eY))
      (IsPreStep.comp (pfRootPre P F) hωs (isPreStep_of_isIso (pfRootPre P F) eA))
  haveI hbA : IsIso ((pfRootPre P F).Base eA) :=
    isBaseIsomorphism_of_isIso (pfRootPre P F) eA
  haveI hbiX : IsIso ((pfRootPre P F).Base (inv eX)) :=
    isBaseIsomorphism_of_isIso (pfRootPre P F) (inv eX)
  haveI hbiY : IsIso ((pfRootPre P F).Base (inv eY)) :=
    isBaseIsomorphism_of_isIso (pfRootPre P F) (inv eY)
  haveI hbζA : IsIso ((pfRootPre P F).Base (ζ ≫ eA)) :=
    (IsPreStep.comp (pfRootPre P F) hζs (isPreStep_of_isIso (pfRootPre P F) eA)).2
  haveI hbωA : IsIso ((pfRootPre P F).Base (ω ≫ eA)) :=
    (IsPreStep.comp (pfRootPre P F) hωs (isPreStep_of_isIso (pfRootPre P F) eA)).2
  have hovζ : (Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base (inv eX ≫ ζ ≫ eA))
        hζ'.2) ((pfRootPre P F).Div (inv eX ≫ ζ ≫ eA))
      = (Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base eA) hbA)
        ((Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base ζ) hζs.2)
          ((pfRootPre P F).Div ζ)) :=
    (overVal_iso_comp (pfRootPre P F) (inv eX) (ζ ≫ eA) hbζA hbiX hζ'.2).trans
      (overVal_comp_iso (pfRootPre P F) ζ eA hζs.2 hbA hbζA)
  have hovω : (Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base (inv eY ≫ ω ≫ eA))
        hω'.2) ((pfRootPre P F).Div (inv eY ≫ ω ≫ eA))
      = (Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base eA) hbA)
        ((Φ.pfOn (phiSharp P)).map (@inv _ _ _ _ ((pfRootPre P F).Base ω) hωs.2)
          ((pfRootPre P F).Div ω)) :=
    (overVal_iso_comp (pfRootPre P F) (inv eY) (ω ≫ eA) hbωA hbiY hω'.2).trans
      (overVal_comp_iso (pfRootPre P F) ω eA hωs.2 hbA hbωA)
  have hle' : MLe ((Φ.pfOn (phiSharp P)).map
        (@inv _ _ _ _ ((pfRootPre P F).Base (inv eY ≫ ω ≫ eA)) hω'.2)
        ((pfRootPre P F).Div (inv eY ≫ ω ≫ eA)))
      ((Φ.pfOn (phiSharp P)).map
        (@inv _ _ _ _ ((pfRootPre P F).Base (inv eX ≫ ζ ≫ eA)) hζ'.2)
        ((pfRootPre P F).Div (inv eX ≫ ζ ≫ eA))) := by
    rw [hovζ, hovω]
    exact MLe.map _ hle
  obtain ⟨m', hm's, hm'e⟩ := pfRoot_coaPreOver_full_sameRoot hfi G
    (inv eX ≫ ζ ≫ eA) (inv eY ≫ ω ≫ eA) hζ' hω' hle'
  refine ⟨eX ≫ m' ≫ inv eY, ?_, ?_⟩
  · exact IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) eX)
      (IsPreStep.comp (pfRootPre P F) hm's (isPreStep_of_isIso (pfRootPre P F) (inv eY)))
  · have h3 := congrArg (fun t => eX ≫ t ≫ inv eA) hm'e
    simp only [Category.assoc, IsIso.hom_inv_id, Category.comp_id,
      IsIso.hom_inv_id_assoc] at h3
    simpa only [Category.assoc] using h3

/-- ★★★★**(iii)(d) 後置の充満性**。 -/
theorem pfRoot_coaPreOver_full (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P)
    (A : PfRootObj P F) :
    letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
    (coaPreOverFunctor (pfRootPre P F) A).Full := by
  letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
  refine ⟨fun {Z W} h => ?_⟩
  obtain ⟨m, hms, hme⟩ := pfRoot_coaPreOver_full_hom hfi G Z.hom.hom W.hom.hom
    Z.hom.property.2 W.hom.property.2 (leOfHom h.unop)
  refine ⟨Over.homMk (show (⟨Z.left.obj⟩ : WideSubcategory (coaPreProp (pfRootPre P F)))
      ⟶ ⟨W.left.obj⟩ from ⟨m, pfRoot_isCoAngular hfi m, hms⟩)
    (WideSubcategory.hom_ext _ hme), Subsingleton.elim _ _⟩

/-! ## ★9. ★★★★★**[FrdI] `Proposition 3.2`** —— `𝒞^pf` は Frobenioid

★これで (iii)(d) の 6 条(前置・後置 × 忠実・充満・本質的全射)が揃った。

| 条 | 宣言 |
|---|---|
| 前置・忠実 | `coaPreUnder_faithful`(pre-Frobenioid 一般) |
| 前置・充満 | `pfRoot_coaPreUnder_full` |
| 前置・本質的全射 | `pfRoot_coaPreUnder_essSurj` |
| 後置・忠実 | `coaPreOver_faithful`(pre-Frobenioid 一般) |
| 後置・充満 | ★`pfRoot_coaPreOver_full` |
| 後置・本質的全射 | `pfRoot_coaPreOver_essSurj` | -/

/-- ★★★★★**[FrdI] `Proposition 3.2`** —— `𝒞` が Frobenius-isotropic 型の Frobenioid なら
**`𝒞^pf` も Frobenioid**。

★21 条は `pfRootCore`、(iii)(d) の 6 条は上の表のとおり。 -/
theorem pfRoot_frobenioid (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P) :
    Frobenioid (pfRootPre P F) where
  core := pfRootCore hfi
  coaPreUnderEquiv := by
    letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
    exact fun X => ⟨coaPreUnder_faithful (pfRootPre P F) X, pfRoot_coaPreUnder_full hfi G X,
      pfRoot_coaPreUnder_essSurj hfi G X⟩
  coaPreOverEquiv := by
    letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
    exact fun X => ⟨coaPreOver_faithful (pfRootPre P F) (fun φ h => pfRoot_preStepMono φ h) X,
      pfRoot_coaPreOver_full hfi G X, pfRoot_coaPreOver_essSurj hfi G X⟩

end SameRoot

/-! ## ★出典の紐付け(条つき) -/

def pfRoot_coaPreOver_essSurj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58,
    item := "Proposition 3.2, (iii) — (iii)(d) 後置の本質的全射性",
    sectionId := "frdi-prop-3-2" }

def pfRoot_coaPreUnder_full.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58,
    item := "Proposition 3.2, (iii) — (iii)(d) 前置の充満性",
    sectionId := "frdi-prop-3-2" }

def pfRoot_coaPreOver_full.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58,
    item := "Proposition 3.2, (iii) — (iii)(d) 後置の充満性",
    sectionId := "frdi-prop-3-2" }

def pfRoot_frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58,
    item := "Proposition 3.2, (iii)",
    sectionId := "frdi-prop-3-2" }

end ABC3.Found.FrdI
