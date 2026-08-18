/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Dict

/-!
# [FrdI] `Proposition 3.2, (iii)` の「perfect 型」

原文 (FrdI p.59):
> (iii) The category Cpf, equipped with the functor Cpf → FΦpf of the diagram of

原文 (FrdI p.59):
> (i), is a Frobenioid of perfect and isotropic type. Moreover, there is a natural

`Frobenioid` は `Prop32Equiv.lean` の `pfRoot_frobenioid` で、
`isotropic 型` は `Prop32Frob.lean` の `pfRoot_isOfIsotropicType` で済んでいる。
★このファイルは **perfect 型**(`Definition 1.2, (iv)`)を扱う。

## ★`IsPerfectObj` の 2 条(`MorphismTypes.lean`)

原文 (FrdI p.23):
> C such that for every n ∈N≥1, it holds that every B ∈Ob(C) base-isomorphic to

(a) 各次数 `n` について、`C` と底同型な対象 `B` は**次数 `n` の Frobenius 型射の
    終域として現れる**。
(b) 次数 `n` の Frobenius 型射 `φ₁ : B₁ ⟶ B₁'`・`φ₂ : B₂ ⟶ B₂'` と
    pre-step `ψ' : B₁' ⟶ B₂'` に対し、`φ₁ ≫ ψ' = ψ ≫ φ₂` を満たす
    **pre-step `ψ : B₁ ⟶ B₂` が一意に存在する**。

## ★★★測定(2026-08-18)—— (b) の一意性は `𝒞` では成り立たない

★`𝒞` では Frobenius 型射は **mono とは限らない**。
`χ₁ ≫ φ = χ₂ ≫ φ` からは(次数・底・零因子が一致するので)
**(vi) `faithfulUpToUnits` により単元のぶんだけのずれ**しか出ない。
★捻れ積 `𝔽_Φ ⋉ G` の模型では `u^n = 1` なる単元 `u` について
`φ ∘ u = φ` が実際に起き、mono は破れる。

★★★**ところが `𝒞^pf` では破れない。** 添字の遷移が Frobenius 型射なので、
**単元は遷移で `u ↦ u^d` と送られ、捻れのぶんが潰れる**。
これを補題の形にしたのが

  ★`homPf_cancel_frobType_idx` —— **Frobenius 型射を添字の遷移に吸収する**。

`ψ` を第 2 脚に持つ遷移 `u` を作ると `φ ≫ ψ = a ≫ idxTransport u φ` であり、
`a` は(pre-Frobenioid の全射性 `totEpiC` で)epi。よって
`φ ≫ ψ = φ' ≫ ψ` から遷移先で `idxTransport u φ = idxTransport u φ'` になる。

★これが `pfRoot_frobTypeMono`(**`𝒞^pf` の Frobenius 型射は mono**)を与え、
(b) の**一意性**が落ちる。

## ★(a) は「根を `n` 倍する」

  `(A, r·n) --Λ_{r·n}(rtExt A n)--> (A^{(n)}, r·n) ≅ (A, r)`

★`𝒞` の `frobDegSurj` は**始域**が自由だが、`𝒞^pf` では**終域**も自由になる
(`pfRoot_frobDegSurj_cod`)。これが perfection の第 1 の効き目である。

## ★★★(b) の**存在** —— 「Frobenius 型射で右から割る」

★★★**`𝒞` ではできない。** `Definition 1.3, (iv)(a)` の 3 分解は Frobenius を
**左**に出すので、右から割るには**添字の自由度**が要る。手は 3 段:

1. 底同型な射は `𝒞` で「**Frobenius 型 ≫ pre-step**」に分解する
   (`exists_frob_preStep_factor`)。3 分解の pull-back 部分が底同型になり、
   `Proposition 1.4, (iii)` で同型に潰れるからである。
2. ★★**`(γ₀, χ, ξ)` を 3 脚添字の遷移として使う** —— `γ₀`(分解の Frobenius 部分)と
   `χ`(割る側の代表元)はどちらも次数 `n` の Frobenius 型射なので、
   これは**正当な遷移**である。遷移先では `θ ↦ π₀ ≫ ξ`・`χ ↦ ξ` になる
   (どちらも epi 性で決まる)。
3. `compRoot_mk3` で `mk(idx12) π₀ ≫ mk(idx23) ξ = mk(idx13) (π₀ ≫ ξ)`。

★これが `pfRoot_frob_factor` であり、`pfRoot_isOfPerfectType` を与える。

## ★残り

`𝒞^pf ≃ (𝒞^pf)^pf`(`Proposition 3.2, (iii)` の後半)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2 v3 u3

section

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★1. 遷移が保つもの -/

/-- ★遷移は Frobenius 型射を保つ。 -/
theorem idxTransport_isFrobeniusType {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (h : IsFrobeniusType P φ) :
    IsFrobeniusType P (idxTransport P F u φ) :=
  prop_1_10_i_frobType_of P F u.right.property.1 u.right.property.2.1
    (idxTransport_spec (F := F) u φ) h

/-- ★遷移は等長性を保つ。 -/
theorem idxTransport_isIsometric {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (h : IsIsometric P φ) :
    IsIsometric P (idxTransport P F u φ) :=
  prop_1_10_i_isometric_of P u.right.property.1.1.2 u.right.property.2.1.1.2
    (idxTransport_spec (F := F) u φ) h

/-- ★遷移は底同型性を保つ。 -/
theorem idxTransport_isBaseIso {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (h : IsBaseIsomorphism P φ) :
    IsBaseIsomorphism P (idxTransport P F u φ) :=
  prop_1_10_i_baseIso_of P u.right.property.1.2 u.right.property.2.1.2
    (idxTransport_spec (F := F) u φ) h

/-! ## ★2. ★★★Frobenius 型射を添字の遷移に吸収する -/

/-- ★★★**`Hom^pf` では Frobenius 型射は右から消せる**。

★★手は「**Frobenius 型射を添字の遷移に吸収する**」こと ——
`ψ` を第 2 脚に持つ遷移 `u` を作ると `idxTransport u φ` は
`φ ≫ ψ = a ≫ idxTransport u φ` を満たす。`a` は(pre-Frobenioid の全射性で)epi
なので、`φ ≫ ψ = φ' ≫ ψ` から遷移先で 2 本が一致する。

★`𝒞` の中では Frobenius 型射は mono とは限らない(単元のぶんだけずれる)。
**`Hom^pf` の中では、そのずれが遷移で潰れる**——これが perfection の効き目である。 -/
theorem homPf_cancel_frobType_idx {A B : C} (Z : IdxPf P F A B)
    (φ φ' : Z.right.obj.1 ⟶ Z.right.obj.2) {E : C} (ψ : Z.right.obj.2 ⟶ E)
    (hψ : IsFrobeniusType P ψ) (h : φ ≫ ψ = φ' ≫ ψ) :
    HomPf.mk Z φ = HomPf.mk Z φ' := by
  obtain ⟨hz1, hz2, hzd⟩ := Z.hom.property
  obtain ⟨W₁, a, ha, had⟩ := F.frobDegSurj Z.right.obj.1 (P.degFr ψ)
  have hdeg : P.degFr (Z.hom.hom.1 ≫ a) = P.degFr (Z.hom.hom.2 ≫ ψ) := by
    rw [P.degFr_comp, P.degFr_comp, had, hzd]
  set u : Z ⟶ idxMk (P := P) (F := F) (Z.hom.hom.1 ≫ a) (Z.hom.hom.2 ≫ ψ)
      (IsFrobeniusType.comp P F hz1 ha) (IsFrobeniusType.comp P F hz2 hψ) hdeg :=
    Under.homMk (show Z.right ⟶ (⟨(W₁, E)⟩ : BiFr P F) from ⟨(a, ψ), ha, hψ, had⟩)
      (WideSubcategory.hom_ext _ rfl) with hu
  have h1 : φ ≫ ψ = a ≫ idxTransport P F u φ := idxTransport_spec (F := F) u φ
  have h2 : φ' ≫ ψ = a ≫ idxTransport P F u φ' := idxTransport_spec (F := F) u φ'
  haveI : Epi a := P.totEpiC _ _ _
  have hkey : idxTransport P F u φ = idxTransport P F u φ' :=
    (cancel_epi a).mp (h1.symm.trans (h.trans h2))
  rw [← HomPf.mk_map u φ, ← HomPf.mk_map u φ', hkey]

set_option maxHeartbeats 1000000 in
/-- ★★★**`Hom^pf` の中で Frobenius 型射は消せる**(3 脚の添字版)。

★`homPf_cancel_preStep` と同じ骨組みだが、最後の一手が違う ——
`𝒞` では Frobenius 型射は mono とは限らないので、
代わりに**添字を isotropic まで押し上げてから**
`homPf_cancel_frobType_idx` を当てる。 -/
theorem homPf_cancel_frobType (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (V : IdxPf3 P F A B E)
    (φ φ' : V.right.obj.1 ⟶ V.right.obj.2.1) (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2)
    (hψi : IsIsometric P ψ) (hψb : IsBaseIsomorphism P ψ)
    (h : HomPf.mk ((idx13 P F A B E).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F A B E).obj V) (φ' ≫ ψ)) :
    HomPf.mk ((idx12 P F A B E).obj V) φ = HomPf.mk ((idx12 P F A B E).obj V) φ' := by
  obtain ⟨V', t, t', ht⟩ := HomPf.eq_iff.mp h
  rw [idx_hom_ext t' t] at ht
  obtain ⟨V'', ⟨k⟩⟩ := exists_hom_of_final (idx13 P F A B E) V'
  set s : V ⟶ IsFiltered.max V V'' := IsFiltered.leftToMax V V'' with hs
  set r : V'' ⟶ IsFiltered.max V V'' := IsFiltered.rightToMax V V'' with hr
  have hm : t ≫ k ≫ (idx13 P F A B E).map r = (idx13 P F A B E).map s :=
    idx_hom_ext _ _
  have hA : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have hB : idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ' ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have key : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ) := by
    rw [hA, hB, ht]
  have hcan : idxTransport P F ((idx12 P F A B E).map s) φ
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ
      = idxTransport P F ((idx12 P F A B E).map s) φ'
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ := by
    rw [idxTransport_comp_pair (F := F) s φ ψ, idxTransport_comp_pair (F := F) s φ' ψ]
    exact key
  -- ★添字を「第 2 脚が isotropic」まで押し上げる
  obtain ⟨M, w, hiso⟩ := exists_idx3_isotropic2 (F := F) hfi (IsFiltered.max V V'')
  have p1 := idxTransport_comp_pair (F := F) w
    (idxTransport P F ((idx12 P F A B E).map s) φ)
    (idxTransport P F ((idx23 P F A B E).map s) ψ)
  have p2 := idxTransport_comp_pair (F := F) w
    (idxTransport P F ((idx12 P F A B E).map s) φ')
    (idxTransport P F ((idx23 P F A B E).map s) ψ)
  have hcanw := p1.trans ((congrArg (idxTransport P F ((idx13 P F A B E).map w)) hcan).trans
    p2.symm)
  have hfrob : IsFrobeniusType P
      (idxTransport P F ((idx23 P F A B E).map w)
        (idxTransport P F ((idx23 P F A B E).map s) ψ)) :=
    (isFrobeniusType_iff_of_isotropic F hiso _).mpr
      ⟨idxTransport_isIsometric _ _ (idxTransport_isIsometric _ _ hψi),
        idxTransport_isBaseIso _ _ (idxTransport_isBaseIso _ _ hψb)⟩
  have hend := homPf_cancel_frobType_idx ((idx12 P F A B E).obj M) _ _ _ hfrob hcanw
  rw [← HomPf.mk_map ((idx12 P F A B E).map s) φ, ← HomPf.mk_map ((idx12 P F A B E).map s) φ',
    ← HomPf.mk_map ((idx12 P F A B E).map w) (idxTransport P F ((idx12 P F A B E).map s) φ),
    ← HomPf.mk_map ((idx12 P F A B E).map w) (idxTransport P F ((idx12 P F A B E).map s) φ')]
  exact hend

set_option maxHeartbeats 1000000 in
/-- ★★★★**`𝒞^pf` の Frobenius 型射は monomorphism**。

★★これは `𝒞` では**成り立たない**——`𝒞` では単元のぶんだけずれる。
`𝒞^pf` では、そのずれが添字の遷移で潰れる(`homPf_cancel_frobType`)。
★これが `Definition 1.2, (iv)` の (b) の**一意性**を与える。 -/
theorem pfRoot_frobTypeMono (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (f : X ⟶ Y) (hf : IsFrobeniusType (pfRootPre P F) f) : Mono f := by
  refine ⟨fun {W} u v huv => ?_⟩
  obtain ⟨V, φ, φ', ψ, hurep, hvrep, hfrep, hu, hv⟩ := compRoot_rep_pairL (F := F) u v f
  have hψ : IsIsometric P ψ ∧ IsBaseIsomorphism P ψ := by
    refine (isFrobeniusType_mk_iff (X := X) (Y := Y) hfi
      ((pushIdx (F := F)
        (rtLift P F X.obj (show Y.root * W.root = W.root * Y.root from mul_comm _ _))
        (rtLift_frobType P F X.obj _)
        (rtLift P F Y.obj (show X.root * W.root = W.root * X.root from mul_comm _ _))
        (rtLift_frobType P F Y.obj _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx23 P F _ _ _).obj V)) ψ).mp ?_
    rw [← rtRootIso_hom_mk (F := F) X.obj Y.obj _ _ ((idx23 P F _ _ _).obj V) ψ, ← hfrep]
    exact hf
  have h0 : (rtRootIso P F W.obj Y.obj
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
        (show X.root * W.root = X.root * W.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ))
      = (rtRootIso P F W.obj Y.obj
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
        (show X.root * W.root = X.root * W.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ)) := by
    rw [← hu, ← hv]; exact huv
  have h1 : HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ) :=
    ((Iso.hom_inv_id_apply _ _).symm.trans (congrArg _ h0)).trans (Iso.hom_inv_id_apply _ _)
  rw [hurep, hvrep, homPf_cancel_frobType (F := F) hfi V φ φ' ψ hψ.1 hψ.2 h1]

/-! ## ★3. `Definition 1.2, (iv)` の (a) —— 終域としての存在 -/

/-- ★★★**`𝒞^pf` では、どの対象も各次数の Frobenius 型射の「終域」として現れる**。

★★これが perfection の第 1 の効き目 —— `𝒞` では `frobDegSurj`(始域が自由)しか
無いが、`𝒞^pf` では**根を `n` 倍した対象から**同じ次数の Frobenius 型射が来る:

  `(A, r·n) --Λ_{r·n}(rtExt A n)--> (A^{(n)}, r·n) ≅ (A, r)`. -/
theorem pfRoot_frobDegSurj_cod (hfi : IsOfFrobeniusIsotropicType P) (B : PfRootObj P F)
    (n : ℕ+) :
    ∃ (B₀ : PfRootObj P F) (φ : B₀ ⟶ B),
      IsFrobeniusType (pfRootPre P F) φ ∧ (pfRootPre P F).degFr φ = n := by
  obtain ⟨e, he⟩ := pfRoot_exists_iso_root (F := F) B.obj B.root n (B.root * n) rfl
  haveI := he
  refine ⟨⟨B.obj, B.root * n⟩,
    lamHom (F := F) (B.root * n) (rtExt P F B.obj n) ≫ inv e, ?_, ?_⟩
  · exact IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi)
      (lamHom_isFrobeniusType (F := F) hfi (B.root * n) (rtExt P F B.obj n)
        (rtExt_frobType P F B.obj n))
      (isFrobeniusType_of_isIso (pfRootPre P F) (inv e))
  · rw [(pfRootPre P F).degFr_comp, lamHom_degFr, rtExt_degFr,
      show (pfRootPre P F).degFr (inv e) = 1 from isLinear_of_isIso (pfRootPre P F) (inv e),
      one_mul]

/-- ★★**(b) の一意性** —— pre-step は Frobenius 型射との合成で決まる。 -/
theorem pfRoot_preStep_uniq_of_comp_frobType (hfi : IsOfFrobeniusIsotropicType P)
    {B₁ B₂ B₂' : PfRootObj P F} (φ₂ : B₂ ⟶ B₂')
    (hφ₂ : IsFrobeniusType (pfRootPre P F) φ₂) (ψ ψ' : B₁ ⟶ B₂)
    (h : ψ ≫ φ₂ = ψ' ≫ φ₂) : ψ = ψ' :=
  haveI := pfRoot_frobTypeMono (F := F) hfi φ₂ hφ₂
  (cancel_mono φ₂).mp h

/-! ## ★4. `Definition 1.2, (iv)` の (b) の**存在** —— Frobenius 型射で「割る」

★★★**要点は 3 つ**:

1. 底同型な射は `𝒞` の中で「**Frobenius 型 ≫ pre-step**」に分解する
   (`exists_frob_preStep_factor`)。`Definition 1.3, (iv)(a)` の 3 分解で
   pull-back の部分が底同型になり、`Proposition 1.4, (iii)` で同型に潰れる。
2. ★★**`(γ₀, χ, ξ)` を 3 脚添字の遷移として使う** —— `γ₀`(分解の Frobenius 部分)と
   `χ`(割る側の代表元)はどちらも次数 `n` の Frobenius 型射なので、これは
   **正当な遷移**である。遷移先では `θ ↦ π₀ ≫ ξ`・`χ ↦ ξ` になる。
3. `compRoot_mk3` で `mk(idx12) π₀ ≫ mk(idx23) ξ = mk(idx13) (π₀ ≫ ξ)`。

★これで「Frobenius 型射で右から割る」が `𝒞^pf` で常に可能になる。
`𝒞` では**できない**——`arbFactor` は Frobenius を**左**に出すので、
右から割るには添字の自由度が要る。 -/

/-- ★★**底同型な射は「Frobenius 型 ≫ pre-step」に分解する**。

★`Definition 1.3, (iv)(a)` の 3 分解 `γ ≫ β ≫ α` で、
底同型性から `α`(pull-back)は **base-isomorphism** にもなり、
`Proposition 1.4, (iii)` により**同型**になる。よって `β ≫ α` は pre-step。
★次数は `degFr γ = degFr θ` になる(`β`・`α` は linear)。 -/
theorem exists_frob_preStep_factor (Fc : FrobenioidCore P) {X Y : C} (θ : X ⟶ Y)
    (hb : IsBaseIsomorphism P θ) :
    ∃ (Y₀ : C) (γ : X ⟶ Y₀) (π : Y₀ ⟶ Y),
      IsFrobeniusType P γ ∧ P.degFr γ = P.degFr θ ∧ IsPreStep P π ∧ θ = γ ≫ π := by
  obtain ⟨X', Y', γ, β, α, hfac, hγ, hβ, hα⟩ := Fc.arbFactor θ
  obtain ⟨hαlb, hαlin⟩ := Fc.pullBackLB α hα
  have hdα : P.degFr α = 1 := hαlin
  have hdβ : P.degFr β = 1 := hβ.1
  have hdγ : P.degFr γ = P.degFr θ := by
    rw [hfac, P.degFr_comp, P.degFr_comp, hdα, hdβ, one_mul, one_mul]
  haveI hbθ : IsIso (P.Base θ) := hb
  haveI hbγ : IsIso (P.Base γ) := hγ.2
  haveI hbβ : IsIso (P.Base β) := hβ.2
  have hsq : P.Base γ ≫ P.Base β ≫ P.Base α = P.Base θ := by
    rw [← P.Base_comp, ← P.Base_comp, ← hfac]
  haveI hbα : IsIso (P.Base α) := by
    haveI : IsIso (P.Base γ ≫ P.Base β ≫ P.Base α) := by rw [hsq]; infer_instance
    haveI : IsIso (P.Base β ≫ P.Base α) := IsIso.of_isIso_comp_left (P.Base γ) _
    exact IsIso.of_isIso_comp_left (P.Base β) _
  haveI : IsIso α := prop_1_4_iii P Fc α hαlb ⟨hαlin, hbα⟩
  exact ⟨X', γ, β ≫ α, hγ, hdγ,
    IsPreStep.comp P hβ (isPreStep_of_isIso P α), by rw [hfac]⟩

set_option maxHeartbeats 2000000 in
/-- ★★★★**Frobenius 型射で右から「割る」**(代表元が揃っている場合)。 -/
theorem pfRoot_frob_factor_rep {A B E : C} {r n : ℕ+}
    (k : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) (φ₂ : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩)
    (V : IdxPf3 P F (rtObj P F A (r * r)) (rtObj P F B (r * r)) (rtObj P F E (r * r)))
    (θ : V.right.obj.1 ⟶ V.right.obj.2.2) (χ : V.right.obj.2.1 ⟶ V.right.obj.2.2)
    (hkmk : k = (rtRootIso P F A E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx13 P F _ _ _).obj V) θ))
    (hφmk : φ₂ = (rtRootIso P F B E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx23 P F _ _ _).obj V) χ))
    (hθb : IsBaseIsomorphism P θ) (hθd : P.degFr θ = n)
    (hχf : IsFrobeniusType P χ) (hχd : P.degFr χ = n) :
    ∃ ψ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩,
      IsPreStep (pfRootPre P F) ψ ∧ ψ ≫ φ₂ = k := by
  obtain ⟨Y₀, γ₀, π₀, hγ₀, hγ₀d, hπ₀, hθfac⟩ := exists_frob_preStep_factor F θ hθb
  obtain ⟨E', ξ, hξ, hξd⟩ := F.frobDegSurj V.right.obj.2.2 n
  have hd12 : P.degFr γ₀ = P.degFr χ := by rw [hγ₀d, hθd, hχd]
  have hd23 : P.degFr χ = P.degFr ξ := by rw [hχd, hξd]
  obtain ⟨hv1, hv2, hv3, hv12, hv23⟩ := V.hom.property
  have ha12 : P.degFr (V.hom.hom.1 ≫ γ₀) = P.degFr (V.hom.hom.2.1 ≫ χ) := by
    rw [P.degFr_comp, P.degFr_comp, hd12, hv12]
  have ha23 : P.degFr (V.hom.hom.2.1 ≫ χ) = P.degFr (V.hom.hom.2.2 ≫ ξ) := by
    rw [P.degFr_comp, P.degFr_comp, hd23, hv23]
  set V' : IdxPf3 P F (rtObj P F A (r * r)) (rtObj P F B (r * r)) (rtObj P F E (r * r)) :=
    idxMk3 (P := P) (F := F) (V.hom.hom.1 ≫ γ₀) (V.hom.hom.2.1 ≫ χ) (V.hom.hom.2.2 ≫ ξ)
      (IsFrobeniusType.comp P F hv1 hγ₀) (IsFrobeniusType.comp P F hv2 hχf)
      (IsFrobeniusType.comp P F hv3 hξ) ha12 ha23 with hV'
  set w : V ⟶ V' :=
    Under.homMk (show V.right ⟶ (⟨(Y₀, V.right.obj.2.2, E')⟩ : TriFr P F) from
      ⟨(γ₀, χ, ξ), hγ₀, hχf, hξ, hd12, hd23⟩) (WideSubcategory.hom_ext _ rfl) with hw
  haveI hepγ : Epi γ₀ := P.totEpiC _ _ _
  haveI hepχ : Epi χ := P.totEpiC _ _ _
  have ht13 : idxTransport P F ((idx13 P F _ _ _).map w) θ = π₀ ≫ ξ := by
    refine ((cancel_epi γ₀).mp ?_).symm
    show γ₀ ≫ (π₀ ≫ ξ) = γ₀ ≫ _
    rw [← Category.assoc, ← hθfac]
    exact idxTransport_spec (F := F) ((idx13 P F _ _ _).map w) θ
  have ht23 : idxTransport P F ((idx23 P F _ _ _).map w) χ = ξ :=
    ((cancel_epi χ).mp (idxTransport_spec (F := F) ((idx23 P F _ _ _).map w) χ)).symm
  have hk' : k = (rtRootIso P F A E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V') (π₀ ≫ ξ)) :=
    hkmk.trans (congrArg _
      ((HomPf.mk_map ((idx13 P F _ _ _).map w) θ).symm.trans
        (congrArg (HomPf.mk ((idx13 P F _ _ _).obj V')) ht13)))
  have hφ' : φ₂ = (rtRootIso P F B E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
      (HomPf.mk ((idx23 P F _ _ _).obj V') ξ) :=
    hφmk.trans (congrArg _
      ((HomPf.mk_map ((idx23 P F _ _ _).map w) χ).symm.trans
        (congrArg (HomPf.mk ((idx23 P F _ _ _).obj V')) ht23)))
  refine ⟨(rtRootIso P F A B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx12 P F _ _ _).obj V') π₀), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    exact (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F B (show r * r = r * r from rfl))
        (rtLift_frobType P F B _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx12 P F _ _ _).obj V')) π₀).mpr hπ₀
  · rw [hk', hφ']
    exact compRoot_mk3 (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      (Z := (⟨E, r⟩ : PfRootObj P F)) V'.hom.hom.1 V'.hom.hom.2.1 V'.hom.hom.2.2
      V'.hom.property.1 V'.hom.property.2.1 V'.hom.property.2.2.1
      V'.hom.property.2.2.2.1 V'.hom.property.2.2.2.2 π₀ ξ

set_option maxHeartbeats 2000000 in
/-- ★★★★**Frobenius 型射で右から「割る」**(同根の場合)。 -/
theorem pfRoot_frob_factor_sameRoot (hfi : IsOfFrobeniusIsotropicType P) {A B E : C} {r n : ℕ+}
    (k : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) (φ₂ : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩)
    (hkd : (pfRootPre P F).degFr k = n) (hkb : IsBaseIsomorphism (pfRootPre P F) k)
    (hφ₂ : IsFrobeniusType (pfRootPre P F) φ₂)
    (hφ₂d : (pfRootPre P F).degFr φ₂ = n) :
    ∃ ψ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩,
      IsPreStep (pfRootPre P F) ψ ∧ ψ ≫ φ₂ = k := by
  obtain ⟨V, θ, χ, hV1, hV2, hf, hg⟩ := exists_rep3_cospan_isotropic (F := F) hfi
    ((rtRootIso P F A E (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv k)
    ((rtRootIso P F B E (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ₂)
  have hkmk : k = (rtRootIso P F A E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx13 P F _ _ _).obj V) θ) := by
    rw [← hf, Iso.inv_hom_id_apply]
  have hφmk : φ₂ = (rtRootIso P F B E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx23 P F _ _ _).obj V) χ) := by
    rw [← hg, Iso.inv_hom_id_apply]
  have hθb : IsBaseIsomorphism P θ := by
    refine (isBaseIsomorphism_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨E, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx13 P F _ _ _).obj V)) θ).mp ?_
    rw [hkmk, rtRootIso_hom_mk] at hkb
    exact hkb
  have hθd : P.degFr θ = n := by
    refine (degFr_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨E, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx13 P F _ _ _).obj V)) θ n).mp ?_
    rw [hkmk, rtRootIso_hom_mk] at hkd
    exact hkd
  have hχf : IsFrobeniusType P χ := by
    refine (isFrobeniusType_mk_iff_of_isotropic (X := (⟨B, r⟩ : PfRootObj P F))
      (Y := (⟨E, r⟩ : PfRootObj P F)) hfi
      ((pushIdx (F := F) (rtLift P F B (show r * r = r * r from rfl))
        (rtLift_frobType P F B _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx23 P F _ _ _).obj V)) χ hV2).mp ?_
    rw [hφmk, rtRootIso_hom_mk] at hφ₂
    exact hφ₂
  have hχd : P.degFr χ = n := by
    refine (degFr_mk_iff (X := (⟨B, r⟩ : PfRootObj P F)) (Y := (⟨E, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F B (show r * r = r * r from rfl))
        (rtLift_frobType P F B _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx23 P F _ _ _).obj V)) χ n).mp ?_
    rw [hφmk, rtRootIso_hom_mk] at hφ₂d
    exact hφ₂d
  exact pfRoot_frob_factor_rep (F := F) k φ₂ V θ χ hkmk hφmk hθb hθd hχf hχd

set_option maxHeartbeats 2000000 in
/-- ★★★★**`𝒞^pf` では Frobenius 型射で右から「割れる」**。

★★★これが `Definition 1.2, (iv)` の (b) の**存在**である。
根を揃えて同根版に落とし、代表元の側で分解と遷移を使う。 -/
theorem pfRoot_frob_factor (hfi : IsOfFrobeniusIsotropicType P) {X Y Z : PfRootObj P F}
    {n : ℕ+} (k : X ⟶ Z) (φ₂ : Y ⟶ Z)
    (hkd : (pfRootPre P F).degFr k = n) (hkb : IsBaseIsomorphism (pfRootPre P F) k)
    (hφ₂ : IsFrobeniusType (pfRootPre P F) φ₂)
    (hφ₂d : (pfRootPre P F).degFr φ₂ = n) :
    ∃ ψ : X ⟶ Y, IsPreStep (pfRootPre P F) ψ ∧ ψ ≫ φ₂ = k := by
  obtain ⟨eX, hX⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root (Y.root * Z.root)
    (X.root * Y.root * Z.root) (by ac_rfl)
  obtain ⟨eY, hY⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root (X.root * Z.root)
    (X.root * Y.root * Z.root) (by ac_rfl)
  obtain ⟨eZ, hZ⟩ := pfRoot_exists_iso_root (F := F) Z.obj Z.root (X.root * Y.root)
    (X.root * Y.root * Z.root) (by ac_rfl)
  haveI := hX; haveI := hY; haveI := hZ
  have hkd' : (pfRootPre P F).degFr (inv eX ≫ k ≫ eZ) = n := by
    rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp,
      show (pfRootPre P F).degFr (inv eX) = 1 from isLinear_of_isIso (pfRootPre P F) (inv eX),
      show (pfRootPre P F).degFr eZ = 1 from isLinear_of_isIso (pfRootPre P F) eZ,
      hkd, mul_one, one_mul]
  have hkb' : IsBaseIsomorphism (pfRootPre P F) (inv eX ≫ k ≫ eZ) := by
    haveI : IsIso ((pfRootPre P F).Base (inv eX)) :=
      isBaseIsomorphism_of_isIso (pfRootPre P F) (inv eX)
    haveI : IsIso ((pfRootPre P F).Base eZ) := isBaseIsomorphism_of_isIso (pfRootPre P F) eZ
    haveI : IsIso ((pfRootPre P F).Base k) := hkb
    show IsIso ((pfRootPre P F).Base (inv eX ≫ k ≫ eZ))
    rw [(pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp]
    infer_instance
  have hφ₂' : IsFrobeniusType (pfRootPre P F) (inv eY ≫ φ₂ ≫ eZ) :=
    IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi)
      (isFrobeniusType_of_isIso (pfRootPre P F) (inv eY))
      (IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi) hφ₂
        (isFrobeniusType_of_isIso (pfRootPre P F) eZ))
  have hφ₂d' : (pfRootPre P F).degFr (inv eY ≫ φ₂ ≫ eZ) = n := by
    rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp,
      show (pfRootPre P F).degFr (inv eY) = 1 from isLinear_of_isIso (pfRootPre P F) (inv eY),
      show (pfRootPre P F).degFr eZ = 1 from isLinear_of_isIso (pfRootPre P F) eZ,
      hφ₂d, mul_one, one_mul]
  obtain ⟨ψ', hψ's, hψ'e⟩ := pfRoot_frob_factor_sameRoot (F := F) hfi
    (inv eX ≫ k ≫ eZ) (inv eY ≫ φ₂ ≫ eZ) hkd' hkb' hφ₂' hφ₂d'
  refine ⟨eX ≫ ψ' ≫ inv eY, ?_, ?_⟩
  · exact IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) eX)
      (IsPreStep.comp (pfRootPre P F) hψ's (isPreStep_of_isIso (pfRootPre P F) (inv eY)))
  · have h3 := congrArg (fun t => eX ≫ t ≫ inv eZ) hψ'e
    simp only [Category.assoc, IsIso.hom_inv_id, Category.comp_id,
      IsIso.hom_inv_id_assoc] at h3
    simpa only [Category.assoc] using h3

/-! ## ★5. ★★★★★perfect 型 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`𝒞^pf` の対象はすべて perfect**(`Definition 1.2, (iv)`)。

(a) は `pfRoot_frobDegSurj_cod`、(b) の存在は `pfRoot_frob_factor`、
(b) の一意性は `pfRoot_frobTypeMono`。 -/
theorem pfRoot_isPerfectObj (hfi : IsOfFrobeniusIsotropicType P) (X : PfRootObj P F) :
    IsPerfectObj (pfRootPre P F) X := by
  intro n
  refine ⟨fun B _ => pfRoot_frobDegSurj_cod (F := F) hfi B n, ?_⟩
  intro B₁ B₁' B₂ B₂' φ₁ φ₂ hf₁ hd₁ hf₂ hd₂ _ _ ψ' hψ'
  have hkd : (pfRootPre P F).degFr (φ₁ ≫ ψ') = n := by
    rw [(pfRootPre P F).degFr_comp, hd₁, show (pfRootPre P F).degFr ψ' = 1 from hψ'.1, one_mul]
  have hkb : IsBaseIsomorphism (pfRootPre P F) (φ₁ ≫ ψ') := by
    haveI : IsIso ((pfRootPre P F).Base φ₁) := hf₁.2
    haveI : IsIso ((pfRootPre P F).Base ψ') := hψ'.2
    show IsIso ((pfRootPre P F).Base (φ₁ ≫ ψ'))
    rw [(pfRootPre P F).Base_comp]
    infer_instance
  obtain ⟨ψ, hψs, hψe⟩ := pfRoot_frob_factor (F := F) hfi (φ₁ ≫ ψ') φ₂ hkd hkb hf₂ hd₂
  haveI := pfRoot_frobTypeMono (F := F) hfi φ₂ hf₂
  refine ⟨ψ, ⟨hψs, hψe.symm⟩, ?_⟩
  rintro ψ₀ ⟨-, hψ₀e⟩
  exact (cancel_mono φ₂).mp (hψ₀e.symm.trans hψe.symm)

/-- ★★★★★**[FrdI] `Proposition 3.2, (iii)` の「perfect 型」の条**。

原文 (FrdI p.59):
> (i), is a Frobenioid of perfect and isotropic type. Moreover, there is a natural
-/
theorem pfRoot_isOfPerfectType (hfi : IsOfFrobeniusIsotropicType P) :
    IsOfPerfectType (pfRootPre P F) :=
  fun X => pfRoot_isPerfectObj (F := F) hfi X

/-! ## ★6. `𝒞^pf ≃ (𝒞^pf)^pf` への下ごしらえ —— 次数が割り切れれば割れる

★perfect 型の (b) は **pre-step** についての主張だが、圏同値の充満性には
**任意の射**について「Frobenius 型射で右から割る」ことが要る。
★手は同じで、分解の Frobenius 部分の次数を `n` に切り揃えるところだけ足す
(`exists_frob_factor_of_dvd`)。 -/

/-- ★★**次数が割り切れれば「次数 `n` の Frobenius 型 ≫ 残り」に分解できる**。

★`Definition 1.3, (iv)(a)` の 3 分解で Frobenius 部分の次数は `degFr θ` そのもの。
★それを `frobDegSurj` で作った `次数 n` ＋ `次数 m` の合成と
`frobDegUniq` で突き合わせると、`n` のところで切れる。 -/
theorem exists_frob_factor_of_dvd (Fc : FrobenioidCore P) {X Y : C} (θ : X ⟶ Y) (n m : ℕ+)
    (hd : P.degFr θ = n * m) :
    ∃ (Y₀ : C) (γ : X ⟶ Y₀) (π : Y₀ ⟶ Y),
      IsFrobeniusType P γ ∧ P.degFr γ = n ∧ θ = γ ≫ π := by
  obtain ⟨X', Y', γ, β, α, hfac, hγ, hβ, hα⟩ := Fc.arbFactor θ
  obtain ⟨hαlb, hαlin⟩ := Fc.pullBackLB α hα
  have hdγ : P.degFr γ = n * m := by
    rw [← hd, hfac, P.degFr_comp, P.degFr_comp,
      show P.degFr α = 1 from hαlin, show P.degFr β = 1 from hβ.1, one_mul, one_mul]
  obtain ⟨X₁, γ₁, hγ₁, hγ₁d⟩ := Fc.frobDegSurj X n
  obtain ⟨X₂, γ₂, hγ₂, hγ₂d⟩ := Fc.frobDegSurj X₁ m
  have hdeq : P.degFr (γ₁ ≫ γ₂) = P.degFr γ := by
    rw [P.degFr_comp, hγ₁d, hγ₂d, hdγ, mul_comm]
  obtain ⟨e, he, hee⟩ := Fc.frobDegUniq X X₂ X' (γ₁ ≫ γ₂) γ
    (IsFrobeniusType.comp P Fc hγ₁ hγ₂) hγ hdeq
  refine ⟨X₁, γ₁, γ₂ ≫ e ≫ β ≫ α, hγ₁, hγ₁d, ?_⟩
  rw [hfac, ← hee]
  simp

set_option maxHeartbeats 2000000 in
/-- ★★★★**Frobenius 型射で右から「割る」**(分解を与えた形)。

★`π₀` が pre-step である必要は無い。pre-step のときはそれが上へ持ち上がる。 -/
theorem pfRoot_frob_div_rep {A B E : C} {r n : ℕ+}
    (k : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) (φ₂ : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩)
    (V : IdxPf3 P F (rtObj P F A (r * r)) (rtObj P F B (r * r)) (rtObj P F E (r * r)))
    (θ : V.right.obj.1 ⟶ V.right.obj.2.2) (χ : V.right.obj.2.1 ⟶ V.right.obj.2.2)
    (hkmk : k = (rtRootIso P F A E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx13 P F _ _ _).obj V) θ))
    (hφmk : φ₂ = (rtRootIso P F B E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx23 P F _ _ _).obj V) χ))
    {Y₀ : C} (γ₀ : V.right.obj.1 ⟶ Y₀) (π₀ : Y₀ ⟶ V.right.obj.2.2)
    (hγ₀ : IsFrobeniusType P γ₀) (hγ₀d : P.degFr γ₀ = n) (hθfac : θ = γ₀ ≫ π₀)
    (hχf : IsFrobeniusType P χ) (hχd : P.degFr χ = n) :
    ∃ ψ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩,
      ψ ≫ φ₂ = k ∧ (IsPreStep P π₀ → IsPreStep (pfRootPre P F) ψ) := by
  obtain ⟨E', ξ, hξ, hξd⟩ := F.frobDegSurj V.right.obj.2.2 n
  have hd12 : P.degFr γ₀ = P.degFr χ := by rw [hγ₀d, hχd]
  have hd23 : P.degFr χ = P.degFr ξ := by rw [hχd, hξd]
  obtain ⟨hv1, hv2, hv3, hv12, hv23⟩ := V.hom.property
  have ha12 : P.degFr (V.hom.hom.1 ≫ γ₀) = P.degFr (V.hom.hom.2.1 ≫ χ) := by
    rw [P.degFr_comp, P.degFr_comp, hd12, hv12]
  have ha23 : P.degFr (V.hom.hom.2.1 ≫ χ) = P.degFr (V.hom.hom.2.2 ≫ ξ) := by
    rw [P.degFr_comp, P.degFr_comp, hd23, hv23]
  set V' : IdxPf3 P F (rtObj P F A (r * r)) (rtObj P F B (r * r)) (rtObj P F E (r * r)) :=
    idxMk3 (P := P) (F := F) (V.hom.hom.1 ≫ γ₀) (V.hom.hom.2.1 ≫ χ) (V.hom.hom.2.2 ≫ ξ)
      (IsFrobeniusType.comp P F hv1 hγ₀) (IsFrobeniusType.comp P F hv2 hχf)
      (IsFrobeniusType.comp P F hv3 hξ) ha12 ha23 with hV'
  set w : V ⟶ V' :=
    Under.homMk (show V.right ⟶ (⟨(Y₀, V.right.obj.2.2, E')⟩ : TriFr P F) from
      ⟨(γ₀, χ, ξ), hγ₀, hχf, hξ, hd12, hd23⟩) (WideSubcategory.hom_ext _ rfl) with hw
  haveI hepγ : Epi γ₀ := P.totEpiC _ _ _
  haveI hepχ : Epi χ := P.totEpiC _ _ _
  have ht13 : idxTransport P F ((idx13 P F _ _ _).map w) θ = π₀ ≫ ξ := by
    refine ((cancel_epi γ₀).mp ?_).symm
    show γ₀ ≫ (π₀ ≫ ξ) = γ₀ ≫ _
    rw [← Category.assoc, ← hθfac]
    exact idxTransport_spec (F := F) ((idx13 P F _ _ _).map w) θ
  have ht23 : idxTransport P F ((idx23 P F _ _ _).map w) χ = ξ :=
    ((cancel_epi χ).mp (idxTransport_spec (F := F) ((idx23 P F _ _ _).map w) χ)).symm
  have hk' : k = (rtRootIso P F A E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V') (π₀ ≫ ξ)) :=
    hkmk.trans (congrArg _
      ((HomPf.mk_map ((idx13 P F _ _ _).map w) θ).symm.trans
        (congrArg (HomPf.mk ((idx13 P F _ _ _).obj V')) ht13)))
  have hφ' : φ₂ = (rtRootIso P F B E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
      (HomPf.mk ((idx23 P F _ _ _).obj V') ξ) :=
    hφmk.trans (congrArg _
      ((HomPf.mk_map ((idx23 P F _ _ _).map w) χ).symm.trans
        (congrArg (HomPf.mk ((idx23 P F _ _ _).obj V')) ht23)))
  refine ⟨(rtRootIso P F A B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx12 P F _ _ _).obj V') π₀), ?_, ?_⟩
  · rw [hk', hφ']
    exact compRoot_mk3 (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      (Z := (⟨E, r⟩ : PfRootObj P F)) V'.hom.hom.1 V'.hom.hom.2.1 V'.hom.hom.2.2
      V'.hom.property.1 V'.hom.property.2.1 V'.hom.property.2.2.1
      V'.hom.property.2.2.2.1 V'.hom.property.2.2.2.2 π₀ ξ
  · intro hπ₀
    rw [rtRootIso_hom_mk]
    exact (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F B (show r * r = r * r from rfl))
        (rtLift_frobType P F B _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx12 P F _ _ _).obj V')) π₀).mpr hπ₀

set_option maxHeartbeats 2000000 in
/-- ★★★★**Frobenius 型射で右から割る**(同根、次数は割り切れればよい)。 -/
theorem pfRoot_frob_div_sameRoot (hfi : IsOfFrobeniusIsotropicType P) {A B E : C} {r n m : ℕ+}
    (k : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) (φ₂ : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩)
    (hkd : (pfRootPre P F).degFr k = n * m)
    (hφ₂ : IsFrobeniusType (pfRootPre P F) φ₂)
    (hφ₂d : (pfRootPre P F).degFr φ₂ = n) :
    ∃ ψ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩, ψ ≫ φ₂ = k := by
  obtain ⟨V, θ, χ, hV1, hV2, hf, hg⟩ := exists_rep3_cospan_isotropic (F := F) hfi
    ((rtRootIso P F A E (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv k)
    ((rtRootIso P F B E (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ₂)
  have hkmk : k = (rtRootIso P F A E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx13 P F _ _ _).obj V) θ) := by
    rw [← hf, Iso.inv_hom_id_apply]
  have hφmk : φ₂ = (rtRootIso P F B E (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom (HomPf.mk ((idx23 P F _ _ _).obj V) χ) := by
    rw [← hg, Iso.inv_hom_id_apply]
  have hθd : P.degFr θ = n * m := by
    refine (degFr_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨E, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl))
        (rtLift_frobType P F A _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx13 P F _ _ _).obj V)) θ (n * m)).mp ?_
    rw [hkmk, rtRootIso_hom_mk] at hkd
    exact hkd
  have hχf : IsFrobeniusType P χ := by
    refine (isFrobeniusType_mk_iff_of_isotropic (X := (⟨B, r⟩ : PfRootObj P F))
      (Y := (⟨E, r⟩ : PfRootObj P F)) hfi
      ((pushIdx (F := F) (rtLift P F B (show r * r = r * r from rfl))
        (rtLift_frobType P F B _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx23 P F _ _ _).obj V)) χ hV2).mp ?_
    rw [hφmk, rtRootIso_hom_mk] at hφ₂
    exact hφ₂
  have hχd : P.degFr χ = n := by
    refine (degFr_mk_iff (X := (⟨B, r⟩ : PfRootObj P F)) (Y := (⟨E, r⟩ : PfRootObj P F))
      ((pushIdx (F := F) (rtLift P F B (show r * r = r * r from rfl))
        (rtLift_frobType P F B _) (rtLift P F E (show r * r = r * r from rfl))
        (rtLift_frobType P F E _) (by rw [rtLift_degFr, rtLift_degFr])).obj
          ((idx23 P F _ _ _).obj V)) χ n).mp ?_
    rw [hφmk, rtRootIso_hom_mk] at hφ₂d
    exact hφ₂d
  obtain ⟨Y₀, γ₀, π₀, hγ₀, hγ₀d, hθfac⟩ := exists_frob_factor_of_dvd F θ n m hθd
  obtain ⟨ψ, hψe, -⟩ := pfRoot_frob_div_rep (F := F) k φ₂ V θ χ hkmk hφmk γ₀ π₀ hγ₀ hγ₀d
    hθfac hχf hχd
  exact ⟨ψ, hψe⟩

set_option maxHeartbeats 2000000 in
/-- ★★★★★**`𝒞^pf` では Frobenius 型射で右から割れる**(次数が割り切れれば十分)。 -/
theorem pfRoot_frob_div (hfi : IsOfFrobeniusIsotropicType P) {X Y Z : PfRootObj P F}
    {n m : ℕ+} (k : X ⟶ Z) (φ₂ : Y ⟶ Z)
    (hkd : (pfRootPre P F).degFr k = n * m)
    (hφ₂ : IsFrobeniusType (pfRootPre P F) φ₂)
    (hφ₂d : (pfRootPre P F).degFr φ₂ = n) :
    ∃ ψ : X ⟶ Y, ψ ≫ φ₂ = k := by
  obtain ⟨eX, hX⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root (Y.root * Z.root)
    (X.root * Y.root * Z.root) (by ac_rfl)
  obtain ⟨eY, hY⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root (X.root * Z.root)
    (X.root * Y.root * Z.root) (by ac_rfl)
  obtain ⟨eZ, hZ⟩ := pfRoot_exists_iso_root (F := F) Z.obj Z.root (X.root * Y.root)
    (X.root * Y.root * Z.root) (by ac_rfl)
  haveI := hX; haveI := hY; haveI := hZ
  have hkd' : (pfRootPre P F).degFr (inv eX ≫ k ≫ eZ) = n * m := by
    rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp,
      show (pfRootPre P F).degFr (inv eX) = 1 from isLinear_of_isIso (pfRootPre P F) (inv eX),
      show (pfRootPre P F).degFr eZ = 1 from isLinear_of_isIso (pfRootPre P F) eZ,
      hkd, mul_one, one_mul]
  have hφ₂' : IsFrobeniusType (pfRootPre P F) (inv eY ≫ φ₂ ≫ eZ) :=
    IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi)
      (isFrobeniusType_of_isIso (pfRootPre P F) (inv eY))
      (IsFrobeniusType.comp (pfRootPre P F) (pfRootCore hfi) hφ₂
        (isFrobeniusType_of_isIso (pfRootPre P F) eZ))
  have hφ₂d' : (pfRootPre P F).degFr (inv eY ≫ φ₂ ≫ eZ) = n := by
    rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp,
      show (pfRootPre P F).degFr (inv eY) = 1 from isLinear_of_isIso (pfRootPre P F) (inv eY),
      show (pfRootPre P F).degFr eZ = 1 from isLinear_of_isIso (pfRootPre P F) eZ,
      hφ₂d, mul_one, one_mul]
  obtain ⟨ψ', hψ'e⟩ := pfRoot_frob_div_sameRoot (F := F) hfi
    (inv eX ≫ k ≫ eZ) (inv eY ≫ φ₂ ≫ eZ) hkd' hφ₂' hφ₂d'
  refine ⟨eX ≫ ψ' ≫ inv eY, ?_⟩
  have h3 := congrArg (fun t => eX ≫ t ≫ inv eZ) hψ'e
  simp only [Category.assoc, IsIso.hom_inv_id, Category.comp_id,
    IsIso.hom_inv_id_assoc] at h3
  simpa only [Category.assoc] using h3

end

/-! ## ★7. `𝒞 → 𝒞^pf` が圏同値になる条件(pre-Frobenioid 一般)

★★**perfection が冪等である**ことの中身。3 条:

| 圏同値の条 | 要るもの |
|---|---|
| 忠実 | **Frobenius 型射が mono** —— 遷移は第 2 脚を右から合成することだから |
| 充満 | **Frobenius 型射で右から割れる** —— 代表元 `ψ` に対し `a ≫ ψ` を `b` で割る |
| 本質的全射 | **どの対象も各次数の Frobenius 型射の終域** —— `(X,r)` を `(B₀,1)` に戻す |

★3 条とも `𝒞^pf` については既に示してある。 -/

section Idem

variable {D : Type u} [Category.{v} D] {C' : Type u3} [Category.{v3} C']
  {Φ : MonoidOn.{v, u, w} D} {Q : PreFrobenioid C' Φ} {Fc : FrobenioidCore Q}

/-- ★★**Frobenius 型射が mono なら `Hom_𝒞(A,B) → Hom^pf(A,B)` は単射**。 -/
theorem toHomPf_injective
    (hmono : ∀ {X Y : C'} (f : X ⟶ Y), IsFrobeniusType Q f → Mono f)
    {A B : C'} {u v : A ⟶ B} (h : toHomPf (F := Fc) u = toHomPf (F := Fc) v) : u = v := by
  obtain ⟨V, t, t', ht⟩ := HomPf.eq_iff.mp h
  rw [idx_hom_ext t' t] at ht
  have h1 := idxTransport_spec (F := Fc) t u
  have h2 := idxTransport_spec (F := Fc) t v
  rw [ht] at h1
  haveI : Mono t.right.hom.2 := hmono _ t.right.property.2.1
  exact (cancel_mono t.right.hom.2).mp (h1.trans h2.symm)

/-- ★★**Frobenius 型射で右から割れるなら `Hom_𝒞(A,B) → Hom^pf(A,B)` は全射**。

★代表元 `ψ` の添字の脚 `(a, b)` に対し `a ≫ ψ` を `b` で割ると、
その商が遷移で `ψ` に写る(`a` は epi)。 -/
theorem toHomPf_surjective
    (hdiv : ∀ {X Y Z : C'} {n m : ℕ+} (k : X ⟶ Z) (φ : Y ⟶ Z),
      IsFrobeniusType Q φ → Q.degFr k = n * m → Q.degFr φ = n → ∃ f : X ⟶ Y, f ≫ φ = k)
    {A B : C'} (g : HomPf Q Fc A B) : ∃ f : A ⟶ B, toHomPf (F := Fc) f = g := by
  obtain ⟨Z, ψ, rfl⟩ := HomPf.exists_rep g
  obtain ⟨hz1, hz2, hzd⟩ := Z.hom.property
  have hkd : Q.degFr (Z.hom.hom.1 ≫ ψ) = Q.degFr Z.hom.hom.2 * Q.degFr ψ := by
    rw [Q.degFr_comp, hzd, mul_comm]
  obtain ⟨f, hf⟩ := hdiv (Z.hom.hom.1 ≫ ψ) Z.hom.hom.2 hz2 hkd rfl
  refine ⟨f, ?_⟩
  set ht : (idxOne Q Fc A B) ⟶ Z :=
    Under.homMk (show (idxOne Q Fc A B).right ⟶ Z.right from ⟨(Z.hom.hom.1, Z.hom.hom.2),
      hz1, hz2, hzd⟩) (WideSubcategory.hom_ext _
        (Prod.ext (Category.id_comp _) (Category.id_comp _))) with hht
  have hspec := idxTransport_spec (F := Fc) ht f
  haveI : Epi Z.hom.hom.1 := Q.totEpiC _ _ _
  have hkey : idxTransport Q Fc ht f = ψ :=
    ((cancel_epi Z.hom.hom.1).mp (hf.symm.trans hspec)).symm
  show HomPf.mk (idxOne Q Fc A B) f = _
  rw [← HomPf.mk_map ht f, hkey]

/-- ★**`𝒞 → 𝒞^pf` は忠実**(Frobenius 型射が mono なら)。 -/
theorem toPfRoot_faithful
    (hmono : ∀ {X Y : C'} (f : X ⟶ Y), IsFrobeniusType Q f → Mono f) :
    (toPfRoot Q Fc).Faithful := by
  refine ⟨fun {A B} {f g} h => ?_⟩
  haveI := isIso_rtExt_one Q Fc A
  haveI := isIso_rtExt_one Q Fc B
  have h2 := toHomPf_injective (Fc := Fc) hmono
    (show toHomPf (F := Fc) (inv (rtExt Q Fc A 1) ≫ f ≫ rtExt Q Fc B 1)
      = toHomPf (F := Fc) (inv (rtExt Q Fc A 1) ≫ g ≫ rtExt Q Fc B 1) from h)
  have h3 := congrArg (fun t => rtExt Q Fc A 1 ≫ t ≫ inv (rtExt Q Fc B 1)) h2
  simpa using h3

/-- ★**`𝒞 → 𝒞^pf` は充満**(Frobenius 型射で割れるなら)。 -/
theorem toPfRoot_full
    (hdiv : ∀ {X Y Z : C'} {n m : ℕ+} (k : X ⟶ Z) (φ : Y ⟶ Z),
      IsFrobeniusType Q φ → Q.degFr k = n * m → Q.degFr φ = n → ∃ f : X ⟶ Y, f ≫ φ = k) :
    (toPfRoot Q Fc).Full := by
  refine ⟨fun {A B} h => ?_⟩
  haveI := isIso_rtExt_one Q Fc A
  haveI := isIso_rtExt_one Q Fc B
  obtain ⟨u, hu⟩ := toHomPf_surjective (Fc := Fc) hdiv
    (show HomPf Q Fc (rtObj Q Fc A 1) (rtObj Q Fc B 1) from h)
  refine ⟨rtExt Q Fc A 1 ≫ u ≫ inv (rtExt Q Fc B 1), ?_⟩
  show toHomPf (F := Fc) (inv (rtExt Q Fc A 1)
      ≫ (rtExt Q Fc A 1 ≫ u ≫ inv (rtExt Q Fc B 1)) ≫ rtExt Q Fc B 1) = _
  rw [show inv (rtExt Q Fc A 1) ≫ (rtExt Q Fc A 1 ≫ u ≫ inv (rtExt Q Fc B 1))
      ≫ rtExt Q Fc B 1 = u from by simp]
  exact hu

/-- ★★**`𝒞 → 𝒞^pf` は本質的全射**(どの対象も各次数の Frobenius 型射の終域なら)。

★`(X, r)` に対し `r` 次の Frobenius 型射 `B₀ ⟶ X` を取ると
`frobDegUniq` で `X ≅ rtObj B₀ r` になり、
`pfRoot_exists_iso_root` の `(B₀,1) ≅ (rtObj B₀ r, r)` と繋がる。 -/
theorem toPfRoot_essSurj
    (hcod : ∀ (B : C') (n : ℕ+), ∃ (B₀ : C') (φ : B₀ ⟶ B),
      IsFrobeniusType Q φ ∧ Q.degFr φ = n) :
    (toPfRoot Q Fc).EssSurj := by
  refine ⟨fun W => ?_⟩
  obtain ⟨B₀, φ, hφ, hφd⟩ := hcod W.obj W.root
  obtain ⟨β, hβ, hβe⟩ := Fc.frobDegUniq B₀ (rtObj Q Fc B₀ W.root) W.obj
    (rtExt Q Fc B₀ W.root) φ (rtExt_frobType Q Fc B₀ W.root) hφ (by rw [rtExt_degFr, hφd])
  obtain ⟨e, he⟩ := pfRoot_exists_iso_root (P := Q) (F := Fc) B₀ 1 W.root W.root
    (one_mul _).symm
  exact ⟨B₀, ⟨(@asIso _ _ _ _ e he) ≪≫
    (@asIso _ _ _ _ (lamHom (F := Fc) W.root β) (lamHom_isIso (F := Fc) W.root β hβ))⟩⟩

/-- ★★★★**3 条が揃えば `𝒞 → 𝒞^pf` は圏同値**。 -/
theorem toPfRoot_isEquivalence
    (hmono : ∀ {X Y : C'} (f : X ⟶ Y), IsFrobeniusType Q f → Mono f)
    (hdiv : ∀ {X Y Z : C'} {n m : ℕ+} (k : X ⟶ Z) (φ : Y ⟶ Z),
      IsFrobeniusType Q φ → Q.degFr k = n * m → Q.degFr φ = n → ∃ f : X ⟶ Y, f ≫ φ = k)
    (hcod : ∀ (B : C') (n : ℕ+), ∃ (B₀ : C') (φ : B₀ ⟶ B),
      IsFrobeniusType Q φ ∧ Q.degFr φ = n) :
    (toPfRoot Q Fc).IsEquivalence :=
  ⟨toPfRoot_faithful hmono, toPfRoot_full hdiv, toPfRoot_essSurj hcod⟩

end Idem

section

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★★★★**[FrdI] `Proposition 3.2, (iii)` の後半** —— `𝒞^pf ≃ (𝒞^pf)^pf`。

原文 (FrdI p.59):
> equivalence of categories Cpf →∼ (Cpf)pf.
-/
theorem pfRoot_pf_isEquivalence (hfi : IsOfFrobeniusIsotropicType P) :
    (toPfRoot (pfRootPre P F) (pfRootCore hfi)).IsEquivalence :=
  toPfRoot_isEquivalence
    (fun f hf => pfRoot_frobTypeMono (F := F) hfi f hf)
    (fun k φ hφ hk hφd => pfRoot_frob_div (F := F) hfi k φ hk hφ hφd)
    (fun B n => pfRoot_frobDegSurj_cod (F := F) hfi B n)

end

/-! ## ★出典の紐付け(条つき) -/

def pfRoot_frobDegSurj_cod.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59,
    item := "Proposition 3.2, (iii) — perfect 型の (a)",
    sectionId := "frdi-prop-3-2" }

def pfRoot_frobTypeMono.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59,
    item := "Proposition 3.2, (iii) — perfect 型の (b) の一意性",
    sectionId := "frdi-prop-3-2" }

def pfRoot_isOfPerfectType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59,
    item := "Proposition 3.2, (iii) — perfect 型",
    sectionId := "frdi-prop-3-2" }

def pfRoot_pf_isEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59,
    item := "Proposition 3.2, (iii) — 𝒞^pf ≃ (𝒞^pf)^pf",
    sectionId := "frdi-prop-3-2" }

/-- ★★★★★**[FrdI] Proposition 3.2 全体**。(i)(ii)(iii) がすべて実装された。

| 条 | 内容 | 実装 |
|---|---|---|
| (i) | 1-可換図式と `𝒞^pf` の pre-Frobenioid 構造 | `pfRootPre`(`Prop32.lean`) |
| (ii) | 10 種の射の辞書 | `Prop32Dict.lean` |
| (iii) | Frobenioid | `pfRoot_frobenioid`(`Prop32Equiv.lean`) |
| (iii) | isotropic 型 | `pfRoot_isOfIsotropicType`(`Prop32Frob.lean`) |
| (iii) | ★perfect 型 | `pfRoot_isOfPerfectType` |
| (iii) | ★`𝒞^pf ≃ (𝒞^pf)^pf` | `pfRoot_pf_isEquivalence` | -/
def prop_3_2.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58, item := "Proposition 3.2",
    sectionId := "frdi-prop-3-2" }

end ABC3.Found.FrdI
