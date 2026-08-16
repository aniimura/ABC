import ABC3.Found.FrdI.Prop18

/-!
# [FrdI] Proposition 1.9 —— Isotropic Objects and Isometries

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.30–p.34(**400 dpi 目視確認 2026-08-15**)。

原文 (FrdI p.30):
> Proposition 1.9. (Isotropic Objects and Isometries) Let Φ be a divisorial

## ★この命題の規模(測定)

**7条 (i)–(vii)、主張は約 31**。これまでで最大。

★**3段階に分けて進める**:
* **第1段: (i)(iv)(vi)(vii)** —— 既存の道具だけで届く。★`Proposition 1.7, (i)(v)` が直接効く
* **第2段: (ii)(iii)** —— (iii) が `Definition 1.3, (i), (c)`(`plBkEquiv`)の初使用
* **第3段: (v)** —— 部分圏 `𝒞^istr` / isotropification 関手 / 随伴

★このファイルは**第1段**である。`.src` は (i)–(vii) が全部揃ってから付ける。

## ★★(v) は「再検証」ではなく「移送」

原文 (FrdI p.33):
> tion 1.4, (i); assertion (iv) [cf. also Definition 1.3, (vii), (b)], it follows immediately

★原文は `𝒞^istr` が `Definition 1.3` を満たすことを
「**`𝒞` が Frobenioid であることから**」導く。**21 条を書き直すのではない。**
(この点は着手前の測定で一度取り違え、原文の証明を読んで訂正した。)
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ### ★準備 —— 原文が「easily verified」と呼ぶ2つ

原文 (FrdI p.32):
> Proof. Since pull-backs which are base-isomorphisms are easily verified to be
-/

/-- ★**base-isomorphism である pull-back は同型**。

`Proposition 1.4, (ii)` で pull-back = LB-invertible ∧ linear、
そこに base-isomorphism を足すと **LB-invertible な pre-step** になり、
`Proposition 1.4, (iii)` で同型。

★原文は「easily verified [cf. Remark 1.2.1]」と書くが、
我々の道では `Proposition 1.4` の (ii)(iii) を経由する。 -/
theorem isIso_of_isPullBack_of_isBaseIso (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hpb : IsPullBack P φ) (hbi : IsBaseIsomorphism P φ) : IsIso φ :=
  prop_1_4_iii P F φ ((prop_1_4_ii P F φ).mp hpb).1
    ⟨((prop_1_4_ii P F φ).mp hpb).2, hbi⟩

/-- ★**co-angular な射を「base-isomorphism ≫ 何か」に分解すると、後半も co-angular**。

`β₁` の3分解を `γ ≫ β₁` の3分解に埋め込むだけ ——
`γ` が base-isomorphism なので、co-angular の条件に現れる
「`α` か `γ'` が base-isomorphism」が**そのまま保たれる**。

★(i) の一意性で、`β = γ ≫ β₁` の `β₁` の co-angular 性を出すのに要る。
`Proposition 1.7, (v)` の co-angular pre-step 版は `β` が pre-step のときしか使えず、
ここでは `β` は base-isomorphism でしかない。 -/
theorem isCoAngular_of_comp_left {A Y X : C} (γ : A ⟶ Y) (β₁ : Y ⟶ X)
    (hco : IsCoAngular P (γ ≫ β₁)) (hγ : IsBaseIsomorphism P γ) : IsCoAngular P β₁ := by
  rintro X' Y' γ' β' α' rfl hα'lin hβ'i hβ's hor
  refine hco X' Y' (γ ≫ γ') β' α' (by simp) hα'lin hβ'i hβ's ?_
  rcases hor with h | h
  · exact Or.inl h
  · exact Or.inr (isBaseIsomorphism_comp P hγ h)

/-- ★**isotropic 性は同型で移る**。

`A ≅ A'` で `A'` が isotropic なら `A` も isotropic ——
`A` から出る isometric pre-step `f` に対し `inv(e) ≫ f` を `A'` の isotropy に当てる。 -/
theorem isIsotropic_of_iso {A A' : C} (e : A ≅ A') (h : IsIsotropic P A') :
    IsIsotropic P A := by
  intro Dd f hfi hfs
  haveI : IsIso (e.inv ≫ f) :=
    h Dd (e.inv ≫ f) (IsIsometric.comp P (isIsometric_of_isIso P e.inv) hfi)
      (IsPreStep.comp P (isPreStep_of_isIso P e.inv) hfs)
  have : f = e.hom ≫ (e.inv ≫ f) := by rw [← Category.assoc, e.hom_inv_id, Category.id_comp]
  rw [this]
  infer_instance

/-! ## ★(i) —— base-isomorphism の「co-angular ≫ isometric pre-step」分解

原文 (FrdI p.31):
> (i) Any base-isomorphism φ : A →B of C admits a factorization

原文 (FrdI p.31):
> where α is an isometric pre-step, and β is a co-angular base-isomorphism; this

★原文の `φ = α ◦ β` は「先に `β`」なので、Lean では `φ = β ≫ α`。
-/

/-- **(i) 存在** —— `Proposition 1.7, (ii)` の分解と
`Definition 1.3, (v), (b)` の分解を**繋ぐだけ**。

1. `Proposition 1.7, (ii)`: `φ = γ ≫ β₀`(`γ` Frobenius 型、`β₀` pre-step)
2. `Definition 1.3, (v), (b)`: `β₀ = β₁ ≫ α₁`(`β₁` co-angular pre-step、`α₁` isometric pre-step)
3. `Proposition 1.7, (i)`: `γ ≫ β₁` は co-angular(合成で閉じる)かつ base-isomorphism -/
theorem prop_1_9_i_factor (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hbi : IsBaseIsomorphism P φ) :
    ∃ (X : C) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ (IsCoAngular P β ∧ IsBaseIsomorphism P β) ∧
        (IsIsometric P α ∧ IsPreStep P α) := by
  obtain ⟨X, γ, β₀, rfl, hγ, hβ₀⟩ := (prop_1_7_ii_baseIso_factor P F φ).mp hbi
  obtain ⟨Y, β₁, α₁, rfl, hβ₁c, hβ₁s, hα₁i, hα₁s⟩ := F.preStepFactor β₀ hβ₀
  exact ⟨Y, γ ≫ β₁, α₁, by rw [Category.assoc],
    ⟨F.coAngularComp γ β₁ hγ.1.1 hβ₁c, isBaseIsomorphism_comp P hγ.2 hβ₁s.2⟩, hα₁i, hα₁s⟩

/-- **(i) 一意性** —— 同型 `δ` を除いて一意。

原文 (FrdI p.31):
> factorization is unique, up to replacing the pair (α, β) by a pair of the form (α ◦

★段取り:
1. `β`, `β'` をそれぞれ `Proposition 1.7, (ii)` で `γ ≫ β₁`、`γ' ≫ β₁'` に分解
2. 次数から `deg_Fr γ = deg_Fr γ'` なので、`Definition 1.3, (ii)` の**本質的一意性**で
   同型 `ε` が取れて `γ ≫ ε = γ'`
3. `𝒞` が totally epimorphic なので `γ` を**消去**して
   `β₁ ≫ α = (ε ≫ β₁') ≫ α'`
4. これは `Definition 1.3, (v), (b)` の**分解の一意性**そのもの -/
theorem prop_1_9_i_uniq (F : FrobenioidCore P) {A B : C} (X X' : C)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβc : IsCoAngular P β) (hβb : IsBaseIsomorphism P β)
    (hαi : IsIsometric P α) (hαs : IsPreStep P α)
    (hβc' : IsCoAngular P β') (hβb' : IsBaseIsomorphism P β')
    (hαi' : IsIsometric P α') (hαs' : IsPreStep P α') :
    ∃ δ : X ≅ X', α' = δ.inv ≫ α ∧ β' = β ≫ δ.hom := by
  obtain ⟨Y, γ, β₁, hβeq, hγ, hβ₁⟩ := (prop_1_7_ii_baseIso_factor P F β).mp hβb
  obtain ⟨Y', γ', β₁', hβeq', hγ', hβ₁'⟩ := (prop_1_7_ii_baseIso_factor P F β').mp hβb'
  -- 2. 次数が一致する
  have hdeg : P.degFr γ = P.degFr γ' := by
    have h1 : P.degFr (β ≫ α) = P.degFr γ := by
      rw [hβeq, P.degFr_comp, P.degFr_comp, hβ₁.1, hαs.1]
      rw [one_mul, one_mul]
    have h2 : P.degFr (β' ≫ α') = P.degFr γ' := by
      rw [hβeq', P.degFr_comp, P.degFr_comp, hβ₁'.1, hαs'.1]
      rw [one_mul, one_mul]
    rw [← h1, heq, h2]
  obtain ⟨ε, hεiso, hγε⟩ := F.frobDegUniq A Y Y' γ γ' hγ hγ' hdeg
  haveI := hεiso
  -- 3. `γ` を消去
  haveI : Epi γ := P.totEpiC _ _ γ
  have hcancel : (β₁ ≫ α : Y ⟶ B) = (ε ≫ β₁') ≫ α' := by
    refine (cancel_epi γ).mp ?_
    rw [← Category.assoc, ← hβeq, heq, hβeq', ← hγε]
    simp
  -- 4. `Definition 1.3, (v), (b)` の一意性
  obtain ⟨δ, hδ1, hδ2⟩ := F.preStepFactorUniq X X' β₁ α (ε ≫ β₁') α' hcancel
    (isCoAngular_of_comp_left P γ β₁ (hβeq ▸ hβc) hγ.2) hβ₁ hαi hαs
    (F.coAngularComp ε β₁' (isCoAngular_of_isIso P ε)
      (isCoAngular_of_comp_left P γ' β₁' (hβeq' ▸ hβc') hγ'.2))
    (IsPreStep.comp P (isPreStep_of_isIso P ε) hβ₁') hαi' hαs'
  refine ⟨δ, hδ1, ?_⟩
  rw [hβeq', ← hγε, Category.assoc, hδ2, hβeq, Category.assoc]

/-- **(i) の3つの ⟺** その1 —— `φ` が isometric ⟺ `β` が Frobenius 型。 -/
theorem prop_1_9_i_isometric_iff {A B X : C} (β : A ⟶ X) (α : X ⟶ B)
    (hβc : IsCoAngular P β) (hβb : IsBaseIsomorphism P β)
    (hαi : IsIsometric P α) (hαs : IsPreStep P α) :
    IsIsometric P (β ≫ α) ↔ IsFrobeniusType P β := by
  constructor
  · intro h
    exact ⟨⟨hβc, (prop_1_7_v_isometric P β α h).1⟩, hβb⟩
  · intro h
    exact IsIsometric.comp P h.1.2 hαi

/-- **(i) の3つの ⟺** その2 —— `φ` が co-angular ⟺ `α` が同型。

★`⇒` は **co-angularity を `α := 𝟙` に当てる技**(`Proposition 1.7` で作ったもの)。 -/
theorem prop_1_9_i_coAngular_iff (F : FrobenioidCore P) {A B X : C} (β : A ⟶ X) (α : X ⟶ B)
    (hβc : IsCoAngular P β) (hβb : IsBaseIsomorphism P β)
    (hαi : IsIsometric P α) (hαs : IsPreStep P α) :
    IsCoAngular P (β ≫ α) ↔ IsIso α := by
  constructor
  · intro h
    exact isIso_of_isCoAngular_right P β α h hαi hαs
  · intro h
    haveI := h
    exact F.coAngularComp β α hβc (isCoAngular_of_isIso P α)

/-- **(i) の3つの ⟺** その3 —— `φ` が pull-back ⟺ `φ` が同型。

★`⇒` が「base-isomorphism である pull-back は同型」そのもの。 -/
theorem prop_1_9_i_pullBack_iff (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hbi : IsBaseIsomorphism P φ) : IsPullBack P φ ↔ IsIso φ := by
  refine ⟨fun h => isIso_of_isPullBack_of_isBaseIso P F φ h hbi, fun h => ?_⟩
  haveI := h
  exact isPullBack_of_isIso P φ

/-! ## ★(iv) —— co-angular linear 射は isotropic 性を両向きに移す

原文 (FrdI p.31):
> (iv) Let φ : A →B be a co-angular linear morphism [e.g., a pull-back
-/

/-- **(iv)** —— `φ` が co-angular linear なら、`A` が isotropic ⟺ `B` が isotropic。

★`⇐` が中身。`A` の isotropic hull `α : A ⟶ A'` を取り、`B` が isotropic なので
普遍性で `φ = α ≫ β` と分解する。`α` は isometric pre-step、`β` は linear
(`Proposition 1.7, (v)`)なので、**co-angularity を `γ := 𝟙` に当てる技**で `α` が同型。 -/
theorem prop_1_9_iv (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hco : IsCoAngular P φ) (hlin : IsLinear P φ) :
    IsIsotropic P A ↔ IsIsotropic P B := by
  refine ⟨fun hA => F.isotropicClosed φ hA, fun hB => ?_⟩
  obtain ⟨A', α, hαi, hαs, hA'iso, huniv⟩ := F.isotropicHullExists A
  obtain ⟨β, hβ, -⟩ := huniv B hB φ
  have hβlin : IsLinear P β := (prop_1_7_v_linear P α β (hβ ▸ hlin)).2
  haveI : IsIso α := isIso_of_isCoAngular_left P α β (hβ ▸ hco) hαi hαs hβlin
  exact isIsotropic_of_iso P (asIso α) hA'iso

/-! ## ★(vi)(vii) —— isotropic hull の2つの特徴づけ

原文 (FrdI p.32):
> (vi) A morphism of C is an isotropic hull if and only if its codomain is

原文 (FrdI p.32):
> isotropic, and, moreover, it is minimal-coadjoint to the morphisms with isotropic

原文 (FrdI p.32):
> (vii) A morphism A →B of C is an isometric pre-step if and only if the com-
-/

/-- `𝒞` の射のうち**始域が isotropic** なもののクラス。 -/
def isotropicDomainProp : MorphismProperty C := fun X _ _ => IsIsotropic P X

/-- ★**isotropic hull に同型を後置しても isotropic hull**。 -/
theorem isIsotropicHull_comp_iso {A A' B : C} (α : A ⟶ A') (β : A' ⟶ B)
    (hα : IsIsotropicHull P α) [IsIso β] : IsIsotropicHull P (α ≫ β) := by
  obtain ⟨hαi, hαs, hA'iso, huniv⟩ := hα
  refine ⟨IsIsometric.comp P hαi (isIsometric_of_isIso P β),
    IsPreStep.comp P hαs (isPreStep_of_isIso P β),
    isIsotropic_of_iso P (asIso β).symm hA'iso, ?_⟩
  intro Cc hCc γ
  obtain ⟨δ', hδ', huq⟩ := huniv Cc hCc γ
  refine ⟨inv β ≫ δ', ?_, ?_⟩
  · show γ = (α ≫ β) ≫ (inv β ≫ δ')
    rw [Category.assoc, ← Category.assoc β, IsIso.hom_inv_id, Category.id_comp]
    exact hδ'
  · intro δ'' hδ''
    have hd0 : γ = (α ≫ β) ≫ δ'' := hδ''
    have hd : γ = α ≫ (β ≫ δ'') := by rw [hd0, Category.assoc]
    have hq : β ≫ δ'' = δ' := huq (β ≫ δ'') hd
    rw [← hq, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]

/-- **(vi)** —— isotropic hull ⟺ 終域が isotropic かつ
「始域が isotropic な射」に minimal-coadjoint。

★`⇒` は普遍性で作った `δ` が `β` の**両側逆射**になることを、
`𝒞` が totally epimorphic であること(`φ` と `α` が epi)から出す。
★`⇐` は isotropic hull の存在から `φ = α ≫ β` を作り、
minimal-coadjoint で `β` を同型にする。 -/
theorem prop_1_9_vi (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsIsotropicHull P φ ↔
      (IsIsotropic P B ∧ IsMinimalCoadjoint (isotropicDomainProp P) φ) := by
  constructor
  · rintro ⟨hφi, hφs, hBiso, huniv⟩
    refine ⟨hBiso, ?_⟩
    rintro X α β rfl hX
    obtain ⟨δ, hδ, -⟩ := huniv X hX α
    haveI : Epi (α ≫ β) := P.totEpiC _ _ _
    haveI : Epi α := P.totEpiC _ _ _
    have h1 : δ ≫ β = 𝟙 B := by
      refine (cancel_epi (α ≫ β)).mp ?_
      rw [Category.comp_id, ← Category.assoc, ← hδ]
    have h2 : β ≫ δ = 𝟙 X := by
      refine (cancel_epi α).mp ?_
      rw [Category.comp_id, ← Category.assoc, ← hδ]
    exact ⟨δ, h2, h1⟩
  · rintro ⟨hBiso, hmc⟩
    obtain ⟨A', α, hα⟩ := F.isotropicHullExists A
    obtain ⟨β, hβ, -⟩ := hα.2.2.2 B hBiso φ
    haveI : IsIso β := hmc A' α β hβ hα.2.2.1
    rw [hβ]
    exact isIsotropicHull_comp_iso P α β hα

/-- ★原文が (vii) の証明の最後に置く観察 ——
`γ = α ◦ β` で3射のうち2つが isometric pre-step なら、残りもそう。

★★**測定: 原文の「any two of the three」は必要より強い。**
`γ = β ≫ α` で **`γ` が isometric pre-step なら、`Proposition 1.7, (v)` により
`α` と `β` の両方が同時に出る**。したがって非自明なのは
「`α`, `β` から `γ`」(= 合成で閉じること)だけである。 -/
theorem isometricPreStep_of_comp {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (hi : IsIsometric P (β ≫ α)) (hs : IsPreStep P (β ≫ α)) :
    (IsIsometric P β ∧ IsPreStep P β) ∧ (IsIsometric P α ∧ IsPreStep P α) :=
  ⟨⟨(prop_1_7_v_isometric P β α hi).1, (prop_1_7_v_preStep P β α hs).1⟩,
   (prop_1_7_v_isometric P β α hi).2, (prop_1_7_v_preStep P β α hs).2⟩

/-- **(vii)** —— `φ : A ⟶ B` が isometric pre-step ⟺
isotropic hull `ψ : B ⟶ C` との合成 `φ ≫ ψ` が `A` の isotropic hull。

★`⇐` は**上の観察1本**で終わる(`φ ≫ ψ` が isometric pre-step なら `φ` もそう)。
★`⇒` は `A` の isotropic hull `α` を取り、普遍性で `φ ≫ ψ = α ≫ β` と書いて
`β` が isometric pre-step であることから(`A'` が isotropic なので)`β` を同型にする。 -/
theorem prop_1_9_vii (F : FrobenioidCore P) {A B Cc : C} (φ : A ⟶ B) (ψ : B ⟶ Cc)
    (hψ : IsIsotropicHull P ψ) :
    (IsIsometric P φ ∧ IsPreStep P φ) ↔ IsIsotropicHull P (φ ≫ ψ) := by
  constructor
  · rintro ⟨hφi, hφs⟩
    obtain ⟨A', α, hα⟩ := F.isotropicHullExists A
    obtain ⟨β, hβ, -⟩ := hα.2.2.2 Cc hψ.2.2.1 (φ ≫ ψ)
    have hcomp : IsIsometric P (φ ≫ ψ) ∧ IsPreStep P (φ ≫ ψ) :=
      ⟨IsIsometric.comp P hφi hψ.1, IsPreStep.comp P hφs hψ.2.1⟩
    have hβis := (isometricPreStep_of_comp P α β (hβ ▸ hcomp.1) (hβ ▸ hcomp.2)).2
    haveI : IsIso β := hα.2.2.1 Cc β hβis.1 hβis.2
    rw [hβ]
    exact isIsotropicHull_comp_iso P α β hα
  · intro h
    exact (isometricPreStep_of_comp P φ ψ h.1 h.2.1).1

/-! ## ★第2段 —— (ii)(iii) の対象レベル

原文 (FrdI p.31):
> (ii) Any base-isomorphism φ : A →B of C induces a functor [well-defined

原文 (FrdI p.31):
> (iii) Any pull-back morphism φ : A →B of C induces a functor [well-

★**原文は「関手 [well-defined up to isomorphism]」と言う。**
まず**対象への割り当てとその一意性**を作る —— そこが数学的な中身であり、
`Definition 1.3, (i), (c)` を使うのもそこである。
-/

/-- **(ii) の対象レベル** —— base-isomorphism `φ : A ⟶ B` は、
`A` への isometric pre-step `ε` を `B` への isometric pre-step `α` へ送る。

★中身は **(i) を合成 `ε ≫ φ` に当てるだけ**。`ε` は base-isomorphism、
`φ` も base-isomorphism なので合成も base-isomorphism で、(i) が使える。 -/
theorem prop_1_9_ii_obj (F : FrobenioidCore P) {A B Cc : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    ∃ (X : C) (β : Cc ⟶ X) (α : X ⟶ B),
      ε ≫ φ = β ≫ α ∧ (IsCoAngular P β ∧ IsBaseIsomorphism P β) ∧
        (IsIsometric P α ∧ IsPreStep P α) :=
  prop_1_9_i_factor P F (ε ≫ φ) (isBaseIsomorphism_comp P hεs.2 hφ)

/-- **(iii) の対象レベル** —— ★★**`Definition 1.3, (i), (c)`(`plBkEquiv`)の初使用**。

pull-back `φ : A ⟶ B` と isometric pre-step `δ : D ⟶ B` に対し、
原文が描く可換四角形

```
Cc --γ--> A
|ψ        |φ
D  --δ--> B
```

を作る。★**手順は原文どおり**:

1. `Base(δ)⁻¹ ∘ Base(φ) : A_𝒟 ⟶ D_𝒟` を `𝒟_{D_𝒟}` の対象と見る
2. **`Definition 1.3, (i), (c)` の圏同値の本質的全射性**で、それを実現する
   pull-back `ψ : Cc ⟶ D` を取る
3. **`Definition 1.2, (ii)`(pull-back の定義そのもの)の全単射**で `γ : Cc ⟶ A` を取る
4. `γ` が isometric pre-step であることは `Proposition 1.7, (v)` から出る -/
theorem prop_1_9_iii_lift (F : FrobenioidCore P) {A B Dd : C} (φ : A ⟶ B)
    (hpb : IsPullBack P φ) (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    ∃ (Cc : C) (γ : Cc ⟶ A) (ψ : Cc ⟶ Dd),
      γ ≫ φ = ψ ≫ δ ∧ IsPullBack P ψ ∧ IsIsometric P γ ∧ IsPreStep P γ := by
  haveI := (F.plBkEquiv Dd).essSurj
  haveI hδb : IsIso (P.Base δ) := hδs.2
  obtain ⟨hφlb, hφlin⟩ := (prop_1_4_ii P F φ).mp hpb
  -- 1./2. 圏同値の**本質的全射性**で pull-back `ψ` を取る
  have hiso := (plBkOverFunctor P Dd).objObjPreimageIso
    (Over.mk (P.Base φ ≫ inv (P.Base δ)))
  set Z := (plBkOverFunctor P Dd).objPreimage
    (Over.mk (P.Base φ ≫ inv (P.Base δ))) with hZ
  have hle : IsIso (Over.Hom.left hiso.hom) :=
    inferInstanceAs (IsIso (((Over.forget ((P.toElem.obj Dd).base)).mapIso hiso).hom))
  -- ★`Over` の射の `left` は型が簡約されないので、素の型の射に落としてから使う
  obtain ⟨e, he⟩ : ∃ e : (P.toElem.obj Z.left.obj).base ⟶ (P.toElem.obj A).base,
      e = Over.Hom.left hiso.hom := ⟨Over.Hom.left hiso.hom, rfl⟩
  haveI hei : IsIso e := by rw [he]; exact hle
  have hw : e ≫ (P.Base φ ≫ inv (P.Base δ)) = P.Base Z.hom.hom := by
    rw [he]; exact Over.w hiso.hom
  -- 3. pull-back の定義の全単射
  have hcond : P.Base (Z.hom.hom ≫ δ) = e ≫ P.Base φ := by
    rw [P.Base_comp, ← hw, Category.assoc, Category.assoc, IsIso.inv_hom_id,
      Category.comp_id]
  obtain ⟨γ, hγ⟩ := (hpb Z.left.obj).2 ⟨(Z.hom.hom ≫ δ, e), hcond⟩
  have hp := Subtype.ext_iff.mp hγ
  have h1 : (γ ≫ φ : Z.left.obj ⟶ B) = Z.hom.hom ≫ δ := congrArg Prod.fst hp
  have h2 : P.Base γ = e := congrArg Prod.snd hp
  obtain ⟨hψlb, hψlin⟩ := (prop_1_4_ii P F Z.hom.hom).mp Z.hom.property
  -- 4. `γ` は isometric pre-step
  refine ⟨Z.left.obj, γ, Z.hom.hom, h1, Z.hom.property, ?_, ?_, ?_⟩
  · refine (prop_1_7_v_isometric P γ φ ?_).1
    rw [h1]
    exact IsIsometric.comp P hψlb.2 hδi
  · have hd : P.degFr φ * P.degFr γ = 1 := by
      rw [← P.degFr_comp, h1, P.degFr_comp, hδs.1, hψlin, one_mul]
    exact pnat_right_eq_one hd
  · show IsIso (P.Base γ)
    rw [h2]
    exact hei

/-- **(ii) の射レベル** —— `φ_*` が射に対して何をするか。

`f : C₁ ⟶ C₂`(isometric pre-step、`f ≫ ε₂ = ε₁`)に対し、
`g : X₁ ⟶ X₂`(isometric pre-step、`g ≫ α₂ = α₁`)を作る。

★段取り: `f ≫ β₂` を **(i) でもう一度分解**して `β' ≫ α'` とすると、
`ε₁ ≫ φ = β' ≫ (α' ≫ α₂)` も (i) の形の分解になる。
**(i) の一意性**で同型 `δ` が取れ、`g := δ.hom ≫ α'` が求めるもの。

★`g` の一意性は `imtrPre_hom_uniq`(`α₂` が mono)から出るので、
**これで `φ_*` の射の割り当ては完全に決まる**。 -/
theorem prop_1_9_ii_hom (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    {C₁ C₂ X₁ X₂ : C} (ε₁ : C₁ ⟶ A) (β₁ : C₁ ⟶ X₁) (α₁ : X₁ ⟶ B)
    (ε₂ : C₂ ⟶ A) (β₂ : C₂ ⟶ X₂) (α₂ : X₂ ⟶ B)
    (he₁ : (ε₁ ≫ φ : C₁ ⟶ B) = β₁ ≫ α₁)
    (hβ₁c : IsCoAngular P β₁) (hβ₁b : IsBaseIsomorphism P β₁)
    (hα₁i : IsIsometric P α₁) (hα₁s : IsPreStep P α₁)
    (he₂ : (ε₂ ≫ φ : C₂ ⟶ B) = β₂ ≫ α₂) (hβ₂b : IsBaseIsomorphism P β₂)
    (hα₂i : IsIsometric P α₂) (hα₂s : IsPreStep P α₂)
    (f : C₁ ⟶ C₂) (hfs : IsPreStep P f) (hf : f ≫ ε₂ = ε₁) :
    ∃ g : X₁ ⟶ X₂, (IsIsometric P g ∧ IsPreStep P g) ∧ g ≫ α₂ = α₁ := by
  obtain ⟨Y', β', α', hfac, ⟨hβ'c, hβ'b⟩, hα'i, hα's⟩ :=
    prop_1_9_i_factor P F (f ≫ β₂) (isBaseIsomorphism_comp P hfs.2 hβ₂b)
  have hkey : (β₁ ≫ α₁ : C₁ ⟶ B) = β' ≫ (α' ≫ α₂) := by
    rw [← he₁, ← hf, Category.assoc, he₂, ← Category.assoc, hfac, Category.assoc]
  obtain ⟨δ, hδ1, -⟩ := prop_1_9_i_uniq P F X₁ Y' β₁ α₁ β' (α' ≫ α₂) hkey
    hβ₁c hβ₁b hα₁i hα₁s hβ'c hβ'b
    (IsIsometric.comp P hα'i hα₂i) (IsPreStep.comp P hα's hα₂s)
  refine ⟨δ.hom ≫ α', ⟨IsIsometric.comp P (isIsometric_of_isIso P δ.hom) hα'i,
    IsPreStep.comp P (isPreStep_of_isIso P δ.hom) hα's⟩, ?_⟩
  rw [Category.assoc, hδ1, ← Category.assoc, δ.hom_inv_id, Category.id_comp]

/-- ★(iii) の四角形を2つ与えると、その間の射が**一意に**作れる。

★使うのは **`ψ'` が pull-back であること**(射を作る)と
**`φ` が pull-back であること**(作った射が `γ` と可換であることを保証する)。
★**圏同値は使わない** —— ここは `Definition 1.2, (ii)` の全単射だけである。 -/
theorem prop_1_9_iii_hom (P : PreFrobenioid C Φ) {A B Dd : C} (φ : A ⟶ B)
    (hpb : IsPullBack P φ) (δ : Dd ⟶ B) (hδs : IsPreStep P δ)
    {X Y : C} (a : X ⟶ A) (b : X ⟶ Dd) (a' : Y ⟶ A) (b' : Y ⟶ Dd)
    (ha : a ≫ φ = b ≫ δ) (hb : a' ≫ φ = b' ≫ δ)
    (ha2 : IsPreStep P a') (hb' : IsPullBack P b') :
    ∃ h : X ⟶ Y, h ≫ a' = a ∧ h ≫ b' = b := by
  haveI hia : IsIso (P.Base a') := ha2.2
  haveI hid : IsIso (P.Base δ) := hδs.2
  have hbase : ∀ {Z : C} (u : Z ⟶ A) (v : Z ⟶ Dd), u ≫ φ = v ≫ δ →
      P.Base v = P.Base u ≫ P.Base φ ≫ inv (P.Base δ) := by
    intro Z u v huv
    have := congrArg P.Base huv
    rw [P.Base_comp, P.Base_comp] at this
    rw [← Category.assoc, IsIso.eq_comp_inv]
    exact this.symm
  have hcond : P.Base b = (P.Base a ≫ inv (P.Base a')) ≫ P.Base b' := by
    rw [hbase a b ha, hbase a' b' hb, Category.assoc, ← Category.assoc (inv (P.Base a')),
      IsIso.inv_hom_id, Category.id_comp]
  obtain ⟨h, hh⟩ := (hb' X).2 ⟨(b, P.Base a ≫ inv (P.Base a')), hcond⟩
  have hp := Subtype.ext_iff.mp hh
  have hb1 : (h ≫ b' : X ⟶ Dd) = b := congrArg Prod.fst hp
  have hb2 : P.Base h = P.Base a ≫ inv (P.Base a') := congrArg Prod.snd hp
  refine ⟨h, ?_, hb1⟩
  refine (hpb X).1 (Subtype.ext (Prod.ext ?_ ?_))
  · show (h ≫ a') ≫ φ = a ≫ φ
    rw [Category.assoc, hb, ← Category.assoc, hb1, ← ha]
  · show P.Base (h ≫ a') = P.Base a
    rw [P.Base_comp, hb2, Category.assoc, IsIso.inv_hom_id, Category.comp_id]

/-- **(iii) の一意性** —— 原文の「the **unique** [up to isomorphism] isometric pre-step」。

★`prop_1_9_iii_hom` を両向きに使い、`γ` が pre-step ゆえ **mono**
(`Definition 1.3, (v), (a)`)であることで両側逆射を得る。 -/
theorem prop_1_9_iii_uniq (F : FrobenioidCore P) {A B Dd : C} (φ : A ⟶ B)
    (hpb : IsPullBack P φ) (δ : Dd ⟶ B) (hδs : IsPreStep P δ)
    {C₁ C₂ : C} (γ₁ : C₁ ⟶ A) (ψ₁ : C₁ ⟶ Dd) (γ₂ : C₂ ⟶ A) (ψ₂ : C₂ ⟶ Dd)
    (h₁ : γ₁ ≫ φ = ψ₁ ≫ δ) (hψ₁ : IsPullBack P ψ₁) (hγ₁ : IsPreStep P γ₁)
    (h₂ : γ₂ ≫ φ = ψ₂ ≫ δ) (hψ₂ : IsPullBack P ψ₂) (hγ₂ : IsPreStep P γ₂) :
    ∃ e : C₁ ≅ C₂, e.hom ≫ γ₂ = γ₁ := by
  obtain ⟨h, hh1, hh2⟩ :=
    prop_1_9_iii_hom P φ hpb δ hδs γ₁ ψ₁ γ₂ ψ₂ h₁ h₂ hγ₂ hψ₂
  obtain ⟨h', hh1', hh2'⟩ :=
    prop_1_9_iii_hom P φ hpb δ hδs γ₂ ψ₂ γ₁ ψ₁ h₂ h₁ hγ₁ hψ₁
  haveI : Mono γ₁ := F.preStepMono γ₁ hγ₁
  haveI : Mono γ₂ := F.preStepMono γ₂ hγ₂
  have e1 : h ≫ h' = 𝟙 C₁ := by
    refine (cancel_mono γ₁).mp ?_
    rw [Category.assoc, hh1', hh1, Category.id_comp]
  have e2 : h' ≫ h = 𝟙 C₂ := by
    refine (cancel_mono γ₂).mp ?_
    rw [Category.assoc, hh1, hh1', Category.id_comp]
  exact ⟨⟨h, h', e1, e2⟩, hh1⟩

/-- ★★**(ii)(iii) の「関手」の射の割り当ては一意に決まる**。

`𝒞^imtr-pre_A` の対象は `A` への isometric pre-step `γ` であり、
その間の射 `h` は `h ≫ γ₂ = γ₁` を満たすもの。
★`γ₂` は pre-step ゆえ **mono**(`Definition 1.3, (v), (a)`)なので、
**そのような `h` は高々1つ**である。

★★**測定**: したがって「関手 [well-defined up to isomorphism]」のうち
**選択が要るのは対象の割り当てだけ**で、**射の割り当てと関手則は自動**である。
★これは (v) の isotropification 関手についての見込み(普遍性が `∃!` なので
射の側に選択は要らない)と**同じ形**であり、ここで先に確かめられた。 -/
theorem imtrPre_hom_uniq (F : FrobenioidCore P) {A C₁ C₂ : C}
    (γ₁ : C₁ ⟶ A) (γ₂ : C₂ ⟶ A) (hγ₂ : IsPreStep P γ₂)
    (h h' : C₁ ⟶ C₂) (hh : h ≫ γ₂ = γ₁) (hh' : h' ≫ γ₂ = γ₁) : h = h' := by
  haveI : Mono γ₂ := F.preStepMono γ₂ hγ₂
  exact (cancel_mono γ₂).mp (hh.trans hh'.symm)

/-! ## ★★(ii) の圏同値の要 —— **ここで忠実性を使う**

原文 (FrdI p.32):
> [cf. the second equivalence of categories of Definition 1.3, (iii), (d)]. Thus, by

★原文は `φ_*` が本質的全射であることを言うために、
「同じ不変量 `(Φ(φ))⁻¹(Div(φ))` を持つ2つの co-angular pre-step は
**同型で移り合う**」を `Definition 1.3, (iii), (d)` から取り出す。

★★**それが下の補題であり、`Definition 1.3, (iii), (d)` の忠実性を使う唯一の箇所である。**
充満性で射を両向きに作り、**忠実性**で合成が恒等になることを言う。 -/

section CoaPreIso

variable [MorphismProperty.IsMultiplicative (coaPreProp P)]

/-- ★★**同じ `Div` を持つ2つの co-angular pre-step は同型で移り合う**。

★使う成分は **充満性と忠実性**。★1.8 では (iii)(d) の忠実性を使わなかったが、
**ここでは使う** —— 「1.8 で使わなかっただけ」であって、
`Definition 1.3, (iii), (d)` の忠実性は**不要な条件ではない**。 -/
theorem coaPre_iso_of_div_eq (hequiv : ∀ X : C, (coaPreUnderFunctor P X).IsEquivalence)
    {A : C} (Z W : Under (⟨A⟩ : WideSubcategory (coaPreProp P)))
    (h : P.Div Z.hom.hom = P.Div W.hom.hom) : Nonempty (Z ≅ W) := by
  haveI := (hequiv A).full
  haveI := (hequiv A).faithful
  have hobj : (coaPreUnderFunctor P A).obj Z = (coaPreUnderFunctor P A).obj W :=
    congrArg toOrderCat h
  obtain ⟨f, -⟩ := (coaPreUnderFunctor P A).map_surjective (eqToHom hobj)
  obtain ⟨g, -⟩ := (coaPreUnderFunctor P A).map_surjective (eqToHom hobj.symm)
  haveI := Preorder.subsingleton_hom ((coaPreUnderFunctor P A).obj Z)
    ((coaPreUnderFunctor P A).obj Z)
  haveI := Preorder.subsingleton_hom ((coaPreUnderFunctor P A).obj W)
    ((coaPreUnderFunctor P A).obj W)
  exact ⟨⟨f, g, (coaPreUnderFunctor P A).map_injective (Subsingleton.elim _ _),
    (coaPreUnderFunctor P A).map_injective (Subsingleton.elim _ _)⟩⟩

end CoaPreIso

/-! ## ★`𝒞^imtr-pre` —— isometric pre-step が定める部分圏

原文 (FrdI p.31):
> Write Cimtr-pre ⊆C for the subcategory determined by the isometric pre-steps
-/

instance : (isometricPreStepProp P).ContainsIdentities :=
  ⟨fun A => ⟨isIsometric_id P A, isPreStep_id P A⟩⟩

instance : (isometricPreStepProp P).IsStableUnderComposition :=
  ⟨fun _ _ hf hg => ⟨IsIsometric.comp P hf.1 hg.1, IsPreStep.comp P hf.2 hg.2⟩⟩

instance : (isometricPreStepProp P).IsMultiplicative where

/-- **`𝒞^imtr-pre`** —— isometric pre-step が定める広い部分圏。

★`𝒞^coa-pre`(`Definition 1.3, (iii), (d)`)と同じ形で作れる ——
`ContainsIdentities` と `IsStableUnderComposition` が
`isIsometric_id` / `isPreStep_id` / `IsIsometric.comp` / `IsPreStep.comp` から
**そのまま出る**(co-angular と違って `Definition 1.3` の条項を引かない)。 -/
abbrev ImtrPre : Type u2 := WideSubcategory (isometricPreStepProp P)

/-- ★**`φ_*` の対象への割り当て** —— `Proposition 1.9, (i)` の分解の
「isometric pre-step 側」を取る。

★★**ここで選択が要る。** `prop_1_9_ii_obj` は **`∃` であって `∃!` ではない** ——
終域は `prop_1_9_i_uniq` により**同型を除いてしか決まらない**。
★対照的に**射の割り当てには選択が要らない**(`imtrPre_hom_uniq`)。
`#print axioms` で `Classical.choice` が入るのはこの定義である。 -/
noncomputable def pushObj (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) {Cc : C} (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) : C :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose

/-- ★上で選んだ対象への isometric pre-step。 -/
noncomputable def pushHom (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) {Cc : C} (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    pushObj P F φ hφ ε hεi hεs ⟶ B :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose_spec.choose

/-- ★`pushHom` は本当に isometric pre-step。 -/
theorem pushHom_spec (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) {Cc : C} (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    IsIsometric P (pushHom P F φ hφ ε hεi hεs) ∧
      IsPreStep P (pushHom P F φ hφ ε hεi hεs) :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose_spec.choose_spec.2.2

/-! ## ★第3段 —— (v) の機械: `𝒞^istr` と isotropification

原文 (FrdI p.32):
> (v) Cistr [equipped with the restriction to C of the given functor C →FΦ] is a

## ★★測定: 部分圏の形が第2段と違う

* `𝒞^imtr-pre` は**射**の条件 → `WideSubcategory`(`InducedWideCategory`)
* `𝒞^istr` は**対象**の条件 → `ObjectProperty.FullSubcategory`(`InducedCategory`)

★**第2段で作った部分圏の機械はここでは使えない。** 見込みどおりだった。
★ただし `InducedCategory.homEquiv` があるので、射の扱いは `WideSubcategory` より**軽い**
(性質のフィールドが無い)。

## ★mathlib の実測(S1–S4、2026-08-15)

| 要るもの | mathlib | 判定 |
|---|---|---|
| 対象の条件による充満部分圏 | `ObjectProperty.FullSubcategory` / `.ι` | ★**使う** |
| induced category の射 | `InducedCategory.homEquiv` / `homMk` | ★**使う** |
| hom 同値から左随伴を作る | `Adjunction.leftAdjointOfEquiv` / `adjunctionOfEquivLeft` | ★★**使う** |

★★**`leftAdjointOfEquiv` は関手と随伴を一度に作る。** isotropic hull の普遍性(`∃!`)が
そのまま hom 同値になるので、**関手を手で組み立てる必要がない**。
-/

section Istr

/-- **`𝒞^istr`** —— isotropic な対象の充満部分圏の述語。 -/
def isotropicProp : ObjectProperty C := fun A => IsIsotropic P A

/-- **`𝒞^istr`**。 -/
abbrev Istr : Type u2 := (isotropicProp P).FullSubcategory

variable (F : FrobenioidCore P)

/-- `A` の isotropic hull の終域(選択)。★`isotropicHullExists` は `∃` なので選択が要る。 -/
noncomputable def hullObj (A : C) : C := (F.isotropicHullExists A).choose

/-- `A` から選んだ isotropic hull への射。 -/
noncomputable def hullMap (A : C) : A ⟶ hullObj P F A :=
  (F.isotropicHullExists A).choose_spec.choose

theorem hullMap_spec (A : C) : IsIsotropicHull P (hullMap P F A) :=
  (F.isotropicHullExists A).choose_spec.choose_spec

/-- `𝒞^istr` の対象としての `A^istr`。 -/
noncomputable def hullIstr (A : C) : Istr P :=
  ⟨hullObj P F A, (hullMap_spec P F A).2.2.1⟩

/-- ★★**isotropic hull の普遍性が、そのまま随伴の hom 同値になる**。

`Hom_{𝒞^istr}(A^istr, Y) ≃ Hom_𝒞(A, Y)`(`Y` は isotropic)。
★`∃!` の存在部分が `right_inv`、一意性部分が `left_inv` になる。 -/
noncomputable def hullHomEquiv (A : C) (Y : Istr P) :
    (hullIstr P F A ⟶ Y) ≃ (A ⟶ (isotropicProp P).ι.obj Y) :=
  InducedCategory.homEquiv.trans
    { toFun := fun g => hullMap P F A ≫ g
      invFun := fun h => ((hullMap_spec P F A).2.2.2 Y.obj Y.property h).choose
      left_inv := fun g =>
        (((hullMap_spec P F A).2.2.2 Y.obj Y.property
          (hullMap P F A ≫ g)).choose_spec.2 g rfl).symm
      right_inv := fun h =>
        (((hullMap_spec P F A).2.2.2 Y.obj Y.property h).choose_spec.1).symm }

/-- ★★**isotropification 関手** —— `leftAdjointOfEquiv` が
**hom 同値から関手そのものを作る**。手で組み立てる必要がない。 -/
noncomputable def isotropification : C ⥤ Istr P :=
  Adjunction.leftAdjointOfEquiv (F_obj := hullIstr P F) (G := (isotropicProp P).ι)
    (e := hullHomEquiv P F)
    (fun X _ _ g h => (Category.assoc (hullMap P F X) h.hom g.hom).symm)

/-- ★★**isotropification は包含関手の左随伴**。

原文 (FrdI p.32):
> Bistr forms a left adjoint to the inclusion functor Cistr →C, through which

★**記録(2026-08-15)**: この行の包含記号は、**PDF の描画では `↪`** だが
**`pdftotext` の抽出では `→`** になる。私は一度 PDF 画像で見たとおり `↪` と写して
ゲートに落とされた。★**「PDF 目視」と「抽出テキスト」が食い違う文字がある**という
測定であり、`▷` / `′` / `≠` が**抽出側で拾えない**のとは別の種類である
(あちらは照合不能、こちらは**別の文字に化ける**)。
★照合できる側(抽出)に合わせ、食い違いを事実として書き残す。 -/
noncomputable def isotropificationAdj : isotropification P F ⊣ (isotropicProp P).ι :=
  Adjunction.adjunctionOfEquivLeft _ _

/-- ★**`𝒞^istr` は totally epimorphic** —— 充満部分圏なので `𝒞` からそのまま移る。

★★これが「移送」の**最初の実例**である: `𝒞` の性質を1行で運ぶ。 -/
theorem istr_totEpi : IsTotallyEpimorphic (Istr P) := by
  intro A B f
  refine ⟨fun {Z} g h hgh => ?_⟩
  haveI : Epi f.hom := P.totEpiC _ _ f.hom
  refine InducedCategory.hom_ext ?_
  exact (cancel_epi f.hom).mp (congrArg InducedCategory.Hom.hom hgh)

include F in
/-- ★**`𝒞` が connected なら `𝒞^istr` も connected**（2026-08-16 に追加）。

★`isotropification` が `𝒞` の zigzag を `𝒞^istr` に送り、
どの isotropic 対象 `Y` も `hullMap` によって `Y^istr` とつながる。
★★**同型である必要はない** —— 辺が 1 本あれば連結性には十分である。 -/
theorem isConnected_istr : IsConnected (Istr P) := by
  haveI := P.connectedC
  obtain ⟨A₀⟩ := (inferInstance : Nonempty C)
  refine IsConnected.of_induct (j₀ := hullIstr P F A₀) ?_
  intro p hp0 hstep Y
  have key : ∀ A : C, hullIstr P F A ∈ p :=
    induct_on_objects (J := C) {A | hullIstr P F A ∈ p} hp0
      (fun {A B} f => hstep ((isotropification P F).map f))
  exact (hstep (⟨hullMap P F Y.obj⟩ : Y ⟶ hullIstr P F Y.obj)).mpr (key Y.obj)

/-- ★★**`𝒞^istr` の pre-Frobenioid 構造** —— 原文の
「equipped with the **restriction** to `𝒞` of the given functor `𝒞 → 𝔽_Φ`」。

★★**4フィールドのうち3つが `P` のものそのまま**である。
`totEpiC` だけが1行の議論を要した。**これが「移送」の意味である。** -/
def istrPre : PreFrobenioid (Istr P) Φ where
  toElem := (isotropicProp P).ι ⋙ P.toElem
  divisorial := P.divisorial
  totEpiC := istr_totEpi P
  totEpiD := P.totEpiD
  connectedC := isConnected_istr P F
  connectedD := P.connectedD

/-- ★`isotropification` の射を、`𝒞` の素の射として取り出したもの。

★`FullSubcategory` の射は `InducedCategory.Hom` に包まれており、
その型が簡約されないので、**素の型に落としてから使う**(第2段と同じ定型)。 -/
noncomputable def istrMap {A B : C} (f : A ⟶ B) : hullObj P F A ⟶ hullObj P F B :=
  ((isotropification P F).map f).hom

/-- ★★**`𝒞^istr` の射の性質は `𝒞` のそれと一致する**(充満部分圏だから)。

原文 (FrdI p.32):
> Cistr satisfies one of these properties with respect to Cistr if and only if it does with

★原文が「compatible with the inclusion functor」と言うのはこれで、
**`rfl` である**(`istrPre` の `toElem` が `ι ⋙ P.toElem` だから)。 -/
theorem istr_compat_degFr {X Y : Istr P} (g : X ⟶ Y) :
    (istrPre P F).degFr g = P.degFr g.hom := rfl

theorem istr_compat_Base {X Y : Istr P} (g : X ⟶ Y) :
    (istrPre P F).Base g = P.Base g.hom := rfl

theorem istr_compat_Div {X Y : Istr P} (g : X ⟶ Y) :
    (istrPre P F).Div g = P.Div g.hom := rfl

/-- ★★**isotropification の定義四角形**。

`leftAdjointOfEquiv` の作る `map f` は、**`hullMap A ≫ istrMap f = f ≫ hullMap B`
を満たす唯一の射**である。原文の「the induced [i.e., by the definition of an
"isotropic hull"!] morphism `A^istr → B^istr`」がこれ。 -/
theorem isotropification_square {A B : C} (f : A ⟶ B) :
    hullMap P F A ≫ istrMap P F f = f ≫ hullMap P F B := by
  have h := (hullHomEquiv P F A (hullIstr P F B)).apply_symm_apply
    (f ≫ hullHomEquiv P F B (hullIstr P F B) (𝟙 _))
  show hullHomEquiv P F A (hullIstr P F B) ((isotropification P F).map f) = _
  rw [show (isotropification P F).map f
      = (hullHomEquiv P F A (hullIstr P F B)).symm
        (f ≫ hullHomEquiv P F B (hullIstr P F B) (𝟙 _)) from rfl, h]
  show f ≫ (hullMap P F B ≫ 𝟙 _) = f ≫ hullMap P F B
  rw [Category.comp_id]

/-! ### ★保存される 11 クラスのうち、まず 3 つの成分

★`Base` / `Div` / `deg_Fr` の保存は、上の四角形と **`Remark 1.1.1`**
(合成公式)だけから出る。原文が「[cf. Remark 1.1.1]」と書くのはこれ。 -/

/-- ★**Frobenius 次数を保つ**。四角形の `deg_Fr` 成分と、
`hullMap` が pre-step(次数 1)であることから。 -/
theorem isotropification_degFr {A B : C} (f : A ⟶ B) :
    P.degFr (istrMap P F f) = P.degFr f := by
  have h := congrArg P.degFr (isotropification_square P F f)
  rw [P.degFr_comp, P.degFr_comp, (hullMap_spec P F A).2.1.1,
    (hullMap_spec P F B).2.1.1, mul_one, one_mul] at h
  exact h

/-- ★**base-isomorphism を(両向きに)保つ**。 -/
theorem isotropification_baseIso_iff {A B : C} (f : A ⟶ B) :
    IsIso (P.Base (istrMap P F f)) ↔ IsIso (P.Base f) := by
  haveI hA : IsIso (P.Base (hullMap P F A)) := (hullMap_spec P F A).2.1.2
  haveI hB : IsIso (P.Base (hullMap P F B)) := (hullMap_spec P F B).2.1.2
  have h := congrArg P.Base (isotropification_square P F f)
  rw [P.Base_comp, P.Base_comp] at h
  constructor
  · intro hi
    haveI := hi
    have hf : P.Base f = P.Base (hullMap P F A) ≫ P.Base (istrMap P F f)
        ≫ inv (P.Base (hullMap P F B)) := by
      rw [← Category.assoc, h, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
    rw [hf]
    infer_instance
  · intro hi
    haveI := hi
    have hf : P.Base (istrMap P F f) = inv (P.Base (hullMap P F A)) ≫ P.Base f
        ≫ P.Base (hullMap P F B) := by
      rw [← h, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    rw [hf]
    infer_instance

/-- ★**isometry を(両向きに)保つ**。四角形の `Div` 成分で、
`hullMap` が isometric なので両端の寄与が消え、`Φ.map` の単射性が残る。 -/
theorem isotropification_isometric_iff {A B : C} (f : A ⟶ B) :
    IsIsometric P (istrMap P F f) ↔ IsIsometric P f := by
  have h := congrArg P.Div (isotropification_square P F f)
  rw [P.Div_comp, P.Div_comp, (hullMap_spec P F A).1, (hullMap_spec P F B).1,
    (hullMap_spec P F B).2.1.1] at h
  simp only [smul_zero, add_zero, PNat.one_coe, one_smul, map_zero, zero_add] at h
  constructor
  · intro hi
    show P.Div f = 0
    rw [← h, show P.Div (istrMap P F f) = 0 from hi, map_zero]
  · intro hi
    show P.Div (istrMap P F f) = 0
    refine Φ.map_injective (P.Base (hullMap P F A)) ?_
    rw [h, show P.Div f = 0 from hi, map_zero]

/-- ★★**isotropic な対象の isotropic hull は同型**。

`hullMap` は isometric pre-step で、`A` が isotropic ならそれは定義から同型。

★原文の「The restriction of the isotropification functor to `𝒞^istr` is
**isomorphic to the identity functor**」の核であり、**1行**である。 -/
theorem hullMap_isIso (A : C) (hA : IsIsotropic P A) : IsIso (hullMap P F A) :=
  hA _ (hullMap P F A) (hullMap_spec P F A).1 (hullMap_spec P F A).2.1

include F in
/-- ★**isotropic な対象から出る `𝒞` の射はすべて co-angular**。

`Definition 1.3, (vii), (b)` で isotropy が伝わるので、`Proposition 1.4, (i)` が使える。
★これが「`𝒞^istr` の射が `𝒞` の意味でも co-angular」を与える —— (v) の
pull-back の保存で要る（原文の「pull-back morphisms **relative to `𝒞`**」）。 -/
theorem isCoAngular_of_isotropic_dom {A B : C} (hA : IsIsotropic P A) (f : A ⟶ B) :
    IsCoAngular P f :=
  prop_1_4_i P f (fun X g => F.isotropicClosed g hA)

/-- ★**pre-step を(両向きに)保つ** —— 次数と底の同型性の合わせ技。 -/
theorem isotropification_preStep_iff {A B : C} (f : A ⟶ B) :
    IsPreStep P (istrMap P F f) ↔ IsPreStep P f := by
  constructor
  · intro h
    refine ⟨?_, (isotropification_baseIso_iff P F f).mp h.2⟩
    show P.degFr f = 1
    rw [← isotropification_degFr P F f]
    exact h.1
  · intro h
    refine ⟨?_, (isotropification_baseIso_iff P F f).mpr h.2⟩
    show P.degFr (istrMap P F f) = 1
    rw [isotropification_degFr P F f]
    exact h.1

/-- ★**pull-back を保つ**（`𝒞` の意味で）。

原文 (FrdI p.33):
> pull-back morphisms to morphisms which are pull-back morphisms relative to C,

★`Proposition 1.4, (ii)` で pull-back = co-angular ∧ isometric ∧ linear に分解し、
* co-angular は `A^istr` が isotropic だから自動（上の補題）
* isometric は保存（`isotropification_isometric_iff`）
* linear は次数の保存

の3つを合わせる。★**`Proposition 1.7` の合成補題は要らない。** -/
theorem isotropification_pullBack {A B : C} (f : A ⟶ B) (h : IsPullBack P f) :
    IsPullBack P (istrMap P F f) := by
  obtain ⟨hlb, hlin⟩ := (prop_1_4_ii P F f).mp h
  refine (prop_1_4_ii P F _).mpr ⟨⟨?_, ?_⟩, ?_⟩
  · exact isCoAngular_of_isotropic_dom P F (hullMap_spec P F A).2.2.1 _
  · exact (isotropification_isometric_iff P F f).mpr hlb.2
  · show P.degFr (istrMap P F f) = 1
    rw [isotropification_degFr P F f]
    exact hlin

/-- ★**Frobenius 型を保つ** —— co-angular（自動）＋ isometric ＋ base-isomorphism。 -/
theorem isotropification_frobType {A B : C} (f : A ⟶ B) (h : IsFrobeniusType P f) :
    IsFrobeniusType P (istrMap P F f) :=
  ⟨⟨isCoAngular_of_isotropic_dom P F (hullMap_spec P F A).2.2.1 _,
    (isotropification_isometric_iff P F f).mpr h.1.2⟩,
   (isotropification_baseIso_iff P F f).mpr h.2⟩

theorem istr_isotropic (X : Istr P) : IsIsotropic (istrPre P F) X := by
  intro Dd φ hi hs
  haveI : IsIso φ.hom := X.property Dd.obj φ.hom hi hs
  exact ⟨InducedCategory.homMk (inv φ.hom), InducedCategory.hom_ext (by simp),
    InducedCategory.hom_ext (by simp)⟩

/-- ★★**`𝒞^istr` のすべての射は co-angular**(`Proposition 1.4, (i)`)。

★原文が保存リストで「co-angular morphisms [cf. Proposition 1.4, (i)]」と
括弧書きするのはこれ —— **`𝒞^istr` では co-angular 性が自明になる**ので、
「保存する」は言うまでもない。 -/
theorem istr_coAngular {X Y : Istr P} (g : X ⟶ Y) : IsCoAngular (istrPre P F) g :=
  prop_1_4_i (istrPre P F) g (fun Z _ => istr_isotropic P F Z)

/-! ### ★★21 条の「移送」—— まず辞書と、移送しやすい条から

★`Definition 1.3` の各条を `istrPre P F` について示す。原文は
「[from the fact that `𝒞` is a Frobenioid!]」の一言だが、
**Lean では `Istr P` と `C` の間の辞書を1本ずつ引く**必要がある。
どこまでが本当に「移送」かを条ごとに測る。 -/

include F in
/-- ★辞書: `Istr P` の Frobenius 型は `C` のそれ。

★co-angular は**両側で自動**なので、実質は isometric と base-isomorphism の移送。 -/
theorem istr_frobType_iff {X Y : Istr P} (g : X ⟶ Y) :
    IsFrobeniusType (istrPre P F) g ↔ IsFrobeniusType P g.hom :=
  ⟨fun h => ⟨⟨isCoAngular_of_isotropic_dom P F X.property g.hom, h.1.2⟩, h.2⟩,
   fun h => ⟨⟨istr_coAngular P F g, h.1.2⟩, h.2⟩⟩

include F in
/-- **(iii)(a)** の移送 —— `𝒞^istr` では co-angular が自動なので**自明**。 -/
theorem istr_coAngularComp {X Y Z : Istr P} (ψ : X ⟶ Y) (φ : Y ⟶ Z) :
    IsCoAngular (istrPre P F) ψ → IsCoAngular (istrPre P F) φ →
      IsCoAngular (istrPre P F) (ψ ≫ φ) :=
  fun _ _ => istr_coAngular P F _

include F in
/-- **(iii)(b)** の移送 —— 同上。 -/
theorem istr_coAngularOfPreStep {X Y : Istr P} (α : X ⟶ Y) :
    IsCoAngular (istrPre P F) α → IsPreStep (istrPre P F) α →
      ∀ φ : X ⟶ Y, IsCoAngular (istrPre P F) φ :=
  fun _ _ φ => istr_coAngular P F φ

include F in
/-- **(vii)(b)** の移送 —— `𝒞^istr` の対象はすべて isotropic なので**自明**。 -/
theorem istr_isotropicClosed {X Y : Istr P} (_φ : X ⟶ Y) :
    IsIsotropic (istrPre P F) X → IsIsotropic (istrPre P F) Y :=
  fun _ => istr_isotropic P F Y

include F in
/-- **(vii)(a)** の移送 —— `X` 自身が isotropic なので `𝟙_X` が isotropic hull。 -/
theorem istr_isotropicHullExists (X : Istr P) :
    ∃ (Y : Istr P) (φ : X ⟶ Y), IsIsotropicHull (istrPre P F) φ :=
  ⟨X, 𝟙 X, (istrPre P F).Div_id X, isPreStep_id _ X, istr_isotropic P F X,
    fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm, fun β hβ => by
      have hg : γ = β := by simpa using hβ
      exact hg.symm⟩⟩

include F in
/-- **(v)(a)** の移送 —— `C` の mono 性が充満部分圏へそのまま降りる。 -/
theorem istr_preStepMono {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P F) φ) :
    Mono φ := by
  haveI : Mono φ.hom := F.preStepMono φ.hom hφ
  refine ⟨fun {Z} g h hgh => ?_⟩
  refine InducedCategory.hom_ext ?_
  exact (cancel_mono φ.hom).mp (congrArg InducedCategory.Hom.hom hgh)

include F in
/-- **(ii)** の本質的一意性の移送 —— `C` で得た同型を充満部分圏へ持ち上げるだけ。 -/
theorem istr_frobDegUniq (X Y Z : Istr P) (φ : X ⟶ Y) (ψ : X ⟶ Z)
    (hφ : IsFrobeniusType (istrPre P F) φ) (hψ : IsFrobeniusType (istrPre P F) ψ)
    (hd : (istrPre P F).degFr φ = (istrPre P F).degFr ψ) :
    ∃ β : Y ⟶ Z, IsIso β ∧ φ ≫ β = ψ := by
  obtain ⟨β₀, hβiso, hβ⟩ := F.frobDegUniq X.obj Y.obj Z.obj φ.hom ψ.hom
    ((istr_frobType_iff P F φ).mp hφ) ((istr_frobType_iff P F ψ).mp hψ) hd
  haveI := hβiso
  refine ⟨InducedCategory.homMk β₀, ?_, InducedCategory.hom_ext hβ⟩
  exact ⟨InducedCategory.homMk (inv β₀), InducedCategory.hom_ext (by simp),
    InducedCategory.hom_ext (by simp)⟩

/-- ★★**pull-back の移送(易しい向き)** —— `𝒞` の pull-back は `𝒞^istr` の pull-back。

原文 (FrdI p.33):
> pull-back morphisms relative to C, hence a fortiori, pull-back morphisms relative to Cistr.

★原文の「**a fortiori**」がこれ。`Definition 1.2, (ii)` の全単射は
「すべての `Z ∈ Ob(𝒞)` について」なので、**`Z` を isotropic なものに制限すれば
そのまま成り立つ**。★ただし Lean では `Istr P` と `C` の射の対応
(`InducedCategory.homEquiv`)を挟む必要がある。

★**逆向き**(「`𝒞^istr` の pull-back は `𝒞` の pull-back」)は
`Z` が isotropic でない場合を埋める必要があり、**随伴を使う**。そちらは別に扱う。 -/
theorem istr_isPullBack_of {X Y : Istr P} (g : X ⟶ Y) (h : IsPullBack P g.hom) :
    IsPullBack (istrPre P F) g := by
  intro Z
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have h1 : (f₁ ≫ g : Z ⟶ Y) = f₂ ≫ g := congrArg Prod.fst hp
    have h2 : P.Base f₁.hom = P.Base f₂.hom := congrArg Prod.snd hp
    refine InducedCategory.hom_ext ?_
    refine (h Z.obj).1 (Subtype.ext (Prod.ext ?_ h2))
    exact congrArg InducedCategory.Hom.hom h1
  · rintro ⟨⟨a, b⟩, hab⟩
    obtain ⟨f₀, hf₀⟩ := (h Z.obj).2 ⟨(a.hom, b), hab⟩
    have hp := Subtype.ext_iff.mp hf₀
    have h1 : (f₀ ≫ g.hom : Z.obj ⟶ Y.obj) = a.hom := congrArg Prod.fst hp
    have h2 : P.Base f₀ = b := congrArg Prod.snd hp
    exact ⟨InducedCategory.homMk f₀,
      Subtype.ext (Prod.ext (InducedCategory.hom_ext h1) h2)⟩

include F in
/-- **(iv)(a)** の移送 —— `𝒞` の3分解を `𝒞^istr` へ運ぶ。

★★**中間対象が自動で isotropic になる**のが効いている ——
`Definition 1.3, (vii), (b)`(`isotropicClosed`)により、
**isotropic な対象から出る射の終域はすべて isotropic** なので、
`𝒞` の分解に現れる対象がそのまま `𝒞^istr` の対象になる。 -/
theorem istr_arbFactor {X Y : Istr P} (φ : X ⟶ Y) :
    ∃ (Z W : Istr P) (γ : X ⟶ Z) (β : Z ⟶ W) (α : W ⟶ Y),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (istrPre P F) γ ∧
        IsPreStep (istrPre P F) β ∧ IsPullBack (istrPre P F) α := by
  obtain ⟨Z₀, W₀, γ₀, β₀, α₀, heq, hγ, hβ, hα⟩ := F.arbFactor φ.hom
  have hZ : IsIsotropic P Z₀ := F.isotropicClosed γ₀ X.property
  have hW : IsIsotropic P W₀ := F.isotropicClosed β₀ hZ
  refine ⟨⟨Z₀, hZ⟩, ⟨W₀, hW⟩, InducedCategory.homMk γ₀, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, ?_, hβ,
    istr_isPullBack_of P F _ hα⟩
  exact (istr_frobType_iff P F (X := X) (Y := ⟨Z₀, hZ⟩)
    (InducedCategory.homMk γ₀)).mpr hγ

include F in
/-- **(ii)** の移送 —— 各次数の Frobenius 型射。中間対象は自動で isotropic。 -/
theorem istr_frobDegSurj (X : Istr P) (n : ℕ+) :
    ∃ (Y : Istr P) (φ : X ⟶ Y), IsFrobeniusType (istrPre P F) φ ∧
      (istrPre P F).degFr φ = n := by
  obtain ⟨B₀, φ₀, hφ, hd⟩ := F.frobDegSurj X.obj n
  exact ⟨⟨B₀, F.isotropicClosed φ₀ X.property⟩, InducedCategory.homMk φ₀,
    (istr_frobType_iff P F _).mpr hφ, hd⟩

include F in
/-- **(v)(b)** の移送。 -/
theorem istr_preStepFactor {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P F) φ) :
    ∃ (Z : Istr P) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsCoAngular (istrPre P F) β ∧ IsPreStep (istrPre P F) β ∧
        IsIsometric (istrPre P F) α ∧ IsPreStep (istrPre P F) α := by
  obtain ⟨Z₀, β₀, α₀, heq, hβc, hβs, hαi, hαs⟩ := F.preStepFactor φ.hom hφ
  refine ⟨⟨Z₀, F.isotropicClosed β₀ X.property⟩, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, ?_, hβs, hαi, hαs⟩
  exact istr_coAngular P F _

include F in
/-- **(v)(c)** の移送。 -/
theorem istr_preStepFactor' {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P F) φ) :
    ∃ (Z : Istr P) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsIsometric (istrPre P F) β ∧ IsPreStep (istrPre P F) β ∧
        IsCoAngular (istrPre P F) α ∧ IsPreStep (istrPre P F) α := by
  obtain ⟨Z₀, β₀, α₀, heq, hβi, hβs, hαc, hαs⟩ := F.preStepFactor' φ.hom hφ
  refine ⟨⟨Z₀, F.isotropicClosed β₀ X.property⟩, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, hβi, hβs, ?_, hαs⟩
  exact istr_coAngular P F _

include F in
/-- ★★**isotropic な対象への pull-back の始域は isotropic**。

`Proposition 1.4, (ii)` で pull-back は **co-angular linear**、
そこに **`Proposition 1.9, (iv)`** を当てるだけ。

★★**原文が (v) の最後で「in light of Proposition 1.4, (i); assertion (iv)」と
書く理由がこれである** —— これがあるから `(𝒞^pl-bk)_{A}` の対象が
**そのまま** `((𝒞^istr)^pl-bk)_{A}` の対象になり、
`Definition 1.3, (i), (c)` の圏同値が `𝒞^istr` へ移送できる。

★**3 行**である。 -/
theorem isotropic_dom_of_pullBack {X A : C} (p : X ⟶ A) (hp : IsPullBack P p)
    (hA : IsIsotropic P A) : IsIsotropic P X := by
  obtain ⟨hlb, hlin⟩ := (prop_1_4_ii P F p).mp hp
  exact (prop_1_9_iv P F p hlb.1 hlin).mpr hA

include F in
/-- ★★**pull-back の移送(難しい向き)** —— `𝒞^istr` の pull-back は `𝒞` の pull-back。

`Definition 1.2, (ii)` の全単射は「**すべての `Z ∈ Ob(𝒞)`**」についてだが、
`istrPre` の側は「**isotropic な `Z`**」しか言わない。この差を埋めるのが

★★**isotropification が包含関手の左随伴であること**(`hullHomEquiv`)である。

`Z` の isotropic hull を取れば `Hom_𝒞(Z, X) ≅ Hom_{𝒞^istr}(Z^istr, X)` なので、
一般の `Z` についての全単射性が isotropic な `Z^istr` のそれに帰着する。
`Base` の側は `Base (hullMap Z)` が同型であることで対応する。

★★**原文はこの依存を書いていない**(「it follows immediately」で済ませている)。
★**随伴は (v) の主張の一部であるだけでなく、(v) を証明する道具でもある。** -/
theorem istr_isPullBack_to {X Y : Istr P} (g : X ⟶ Y) (h : IsPullBack (istrPre P F) g) :
    IsPullBack P g.hom := by
  intro Z
  haveI hZb : IsIso (P.Base (hullMap P F Z)) := (hullMap_spec P F Z).2.1.2
  haveI hZe : Epi (hullMap P F Z) := P.totEpiC _ _ _
  set u := hullMap P F Z with hu
  haveI hub : IsIso (P.Base u) := hZb
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have h1 : (f₁ ≫ g.hom : Z ⟶ Y.obj) = f₂ ≫ g.hom := congrArg Prod.fst hp
    have h2 : P.Base f₁ = P.Base f₂ := congrArg Prod.snd hp
    have e₁ : u ≫ ((hullHomEquiv P F Z X).symm f₁).hom = f₁ :=
      (hullHomEquiv P F Z X).apply_symm_apply f₁
    have e₂ : u ≫ ((hullHomEquiv P F Z X).symm f₂).hom = f₂ :=
      (hullHomEquiv P F Z X).apply_symm_apply f₂
    have hgg : (hullHomEquiv P F Z X).symm f₁ = (hullHomEquiv P F Z X).symm f₂ := by
      refine (h (hullIstr P F Z)).1 (Subtype.ext (Prod.ext ?_ ?_))
      · refine InducedCategory.hom_ext ?_
        refine (cancel_epi u).mp ?_
        calc u ≫ (((hullHomEquiv P F Z X).symm f₁) ≫ g).hom
            = (u ≫ ((hullHomEquiv P F Z X).symm f₁).hom) ≫ g.hom :=
              (Category.assoc _ _ _).symm
          _ = f₁ ≫ g.hom := by rw [e₁]
          _ = f₂ ≫ g.hom := h1
          _ = (u ≫ ((hullHomEquiv P F Z X).symm f₂).hom) ≫ g.hom := by rw [e₂]
          _ = u ≫ (((hullHomEquiv P F Z X).symm f₂) ≫ g).hom := Category.assoc _ _ _
      · show P.Base ((hullHomEquiv P F Z X).symm f₁).hom
          = P.Base ((hullHomEquiv P F Z X).symm f₂).hom
        refine (cancel_epi (P.Base u)).mp ?_
        calc P.Base u ≫ P.Base ((hullHomEquiv P F Z X).symm f₁).hom
            = P.Base (u ≫ ((hullHomEquiv P F Z X).symm f₁).hom) := (P.Base_comp _ _).symm
          _ = P.Base f₁ := by rw [e₁]
          _ = P.Base f₂ := h2
          _ = P.Base (u ≫ ((hullHomEquiv P F Z X).symm f₂).hom) := by rw [e₂]
          _ = P.Base u ≫ P.Base ((hullHomEquiv P F Z X).symm f₂).hom := P.Base_comp _ _
    calc f₁ = u ≫ ((hullHomEquiv P F Z X).symm f₁).hom := e₁.symm
      _ = u ≫ ((hullHomEquiv P F Z X).symm f₂).hom := by rw [hgg]
      _ = f₂ := e₂
  · rintro ⟨⟨a, b⟩, hab⟩
    obtain ⟨w, hw1, hw2⟩ := hZb.out
    have ea : u ≫ ((hullHomEquiv P F Z Y).symm a).hom = a :=
      (hullHomEquiv P F Z Y).apply_symm_apply a
    have hcond : (istrPre P F).Base ((hullHomEquiv P F Z Y).symm a)
        = (w ≫ b) ≫ (istrPre P F).Base g := by
      show P.Base ((hullHomEquiv P F Z Y).symm a).hom
        = (w ≫ b) ≫ P.Base g.hom
      have hbase : P.Base u ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom
          = b ≫ P.Base g.hom := by
        rw [← P.Base_comp, ea, hab]
      calc P.Base ((hullHomEquiv P F Z Y).symm a).hom
          = 𝟙 _ ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom := (Category.id_comp _).symm
        _ = (w ≫ P.Base u) ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom :=
              congrArg (fun t => t ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom) hw2.symm
        _ = w ≫ (P.Base u ≫ P.Base ((hullHomEquiv P F Z Y).symm a).hom) :=
              Category.assoc _ _ _
        _ = w ≫ (b ≫ P.Base g.hom) := by rw [hbase]
        _ = (w ≫ b) ≫ P.Base g.hom := (Category.assoc _ _ _).symm
    obtain ⟨f', hf'⟩ := (h (hullIstr P F Z)).2
      ⟨((hullHomEquiv P F Z Y).symm a, w ≫ b), hcond⟩
    have hp := Subtype.ext_iff.mp hf'
    have k1 : (f' ≫ g : hullIstr P F Z ⟶ Y) = (hullHomEquiv P F Z Y).symm a :=
      congrArg Prod.fst hp
    have k2 : P.Base f'.hom = w ≫ b := congrArg Prod.snd hp
    have hk : f'.hom ≫ g.hom = ((hullHomEquiv P F Z Y).symm a).hom :=
      congrArg InducedCategory.Hom.hom k1
    refine ⟨u ≫ f'.hom, Subtype.ext (Prod.ext ?_ ?_)⟩
    · show (u ≫ f'.hom) ≫ g.hom = a
      calc (u ≫ f'.hom) ≫ g.hom = u ≫ (f'.hom ≫ g.hom) := Category.assoc _ _ _
        _ = u ≫ ((hullHomEquiv P F Z Y).symm a).hom := congrArg (fun t => u ≫ t) hk
        _ = a := ea
    · show P.Base (u ≫ f'.hom) = b
      calc P.Base (u ≫ f'.hom) = P.Base u ≫ P.Base f'.hom := P.Base_comp _ _
        _ = P.Base u ≫ (w ≫ b) := congrArg (fun t => P.Base u ≫ t) k2
        _ = (P.Base u ≫ w) ≫ b := (Category.assoc _ _ _).symm
        _ = 𝟙 _ ≫ b := by rw [hw1]
        _ = b := Category.id_comp _

/-! ### ★残りの条の移送 —— 「前向き／後ろ向き」の仕分け

★**前向き**(与えられた対象の間で射を作る・一意性を言う)は、
新しい対象を作らないので**そのまま移送できる**。
★**後ろ向き**(`𝒟` の情報から `𝒞` の対象を作る)は **isotropification が要る**。

| 条 | 向き |
|---|---|
| `pullBackLB` / `preStepFactorUniq` / `preStepFactorUniq'` / `arbFactorUniq` | **前向き** |
| `faithfulUpToUnits` / `otriFwd` / `otriBwd` / `otriBase` | **前向き** |
| `baseSurj` / `preStepSpan` | ★**後ろ向き** |
| `plBkEquiv` | 圏同値(別扱い) |
-/

include F in
/-- **(iv)(b)** の移送 —— ★**逆向きの pull-back 移送が開いた**ので通る。 -/
theorem istr_pullBackLB {X Y : Istr P} (α : X ⟶ Y) (h : IsPullBack (istrPre P F) α) :
    IsLBInvertible (istrPre P F) α ∧ IsLinear (istrPre P F) α := by
  obtain ⟨hlb, hlin⟩ := F.pullBackLB α.hom (istr_isPullBack_to P F α h)
  exact ⟨⟨istr_coAngular P F α, hlb.2⟩, hlin⟩

include F in
/-- **(v)(b)** の一意性の移送 —— `C` で得た同型を充満部分圏へ持ち上げるだけ。 -/
theorem istr_preStepFactorUniq {A B : Istr P} (X X' : Istr P)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβc : IsCoAngular (istrPre P F) β) (hβs : IsPreStep (istrPre P F) β)
    (hαi : IsIsometric (istrPre P F) α) (hαs : IsPreStep (istrPre P F) α)
    (hβc' : IsCoAngular (istrPre P F) β') (hβs' : IsPreStep (istrPre P F) β')
    (hαi' : IsIsometric (istrPre P F) α') (hαs' : IsPreStep (istrPre P F) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨γ₀, h1, h2⟩ := F.preStepFactorUniq X.obj X'.obj β.hom α.hom β'.hom α'.hom
    (congrArg InducedCategory.Hom.hom heq)
    (isCoAngular_of_isotropic_dom P F A.property _) hβs hαi hαs
    (isCoAngular_of_isotropic_dom P F A.property _) hβs' hαi' hαs'
  exact ⟨InducedCategory.isoMk γ₀, InducedCategory.hom_ext h1, InducedCategory.hom_ext h2⟩

include F in
/-- **(v)(c)** の一意性の移送。 -/
theorem istr_preStepFactorUniq' {A B : Istr P} (X X' : Istr P)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβi : IsIsometric (istrPre P F) β) (hβs : IsPreStep (istrPre P F) β)
    (hαc : IsCoAngular (istrPre P F) α) (hαs : IsPreStep (istrPre P F) α)
    (hβi' : IsIsometric (istrPre P F) β') (hβs' : IsPreStep (istrPre P F) β')
    (hαc' : IsCoAngular (istrPre P F) α') (hαs' : IsPreStep (istrPre P F) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨γ₀, h1, h2⟩ := F.preStepFactorUniq' X.obj X'.obj β.hom α.hom β'.hom α'.hom
    (congrArg InducedCategory.Hom.hom heq)
    hβi hβs (isCoAngular_of_isotropic_dom P F X.property _) hαs
    hβi' hβs' (isCoAngular_of_isotropic_dom P F X'.property _) hαs'
  exact ⟨InducedCategory.isoMk γ₀, InducedCategory.hom_ext h1, InducedCategory.hom_ext h2⟩

include F in
/-- **(vi)** の移送 —— `𝒪^×` の元も `C` から持ち上がる。 -/
theorem istr_faithfulUpToUnits {A B : Istr P} (φ ψ : A ⟶ B)
    (hb : BaseEquivalent (istrPre P F) φ ψ) (hm : MetricallyEquivalent (istrPre P F) φ ψ)
    (hφc : IsCoAngular (istrPre P F) φ) (hφs : IsPreStep (istrPre P F) φ)
    (hψc : IsCoAngular (istrPre P F) ψ) (hψs : IsPreStep (istrPre P F) ψ) :
    ∃ α : End B, α ∈ OTimes (istrPre P F) B ∧ φ = ψ ≫ (α : B ⟶ B) := by
  obtain ⟨α₀, hα₀, hφψ⟩ := F.faithfulUpToUnits φ.hom ψ.hom hb hm
    (isCoAngular_of_isotropic_dom P F A.property _) hφs
    (isCoAngular_of_isotropic_dom P F A.property _) hψs
  refine ⟨InducedCategory.homMk α₀, ⟨⟨hα₀.1.1, hα₀.1.2⟩, ?_⟩, InducedCategory.hom_ext hφψ⟩
  obtain ⟨v, hv⟩ := hα₀.2
  haveI : IsIso (α₀ : B.obj ⟶ B.obj) := (isUnit_iff_isIso _).mp hα₀.2
  refine (isUnit_iff_isIso _).mpr ?_
  exact ⟨InducedCategory.homMk (inv (α₀ : B.obj ⟶ B.obj)),
    InducedCategory.hom_ext (by simp), InducedCategory.hom_ext (by simp)⟩

include F in
/-- **(iv)(a)** の一意性の移送 —— 2つの同型をどちらも持ち上げる。 -/
theorem istr_arbFactorUniq {A B : Istr P} (X Y X' Y' : Istr P)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (heq : (γ ≫ β ≫ α : A ⟶ B) = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType (istrPre P F) γ) (hβ : IsPreStep (istrPre P F) β)
    (hα : IsPullBack (istrPre P F) α)
    (hγ' : IsFrobeniusType (istrPre P F) γ') (hβ' : IsPreStep (istrPre P F) β')
    (hα' : IsPullBack (istrPre P F) α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  obtain ⟨δ₀, ε₀, h1, h2, h3⟩ := F.arbFactorUniq X.obj Y.obj X'.obj Y'.obj
    γ.hom β.hom α.hom γ'.hom β'.hom α'.hom
    (congrArg InducedCategory.Hom.hom heq)
    ((istr_frobType_iff P F γ).mp hγ) hβ (istr_isPullBack_to P F α hα)
    ((istr_frobType_iff P F γ').mp hγ') hβ' (istr_isPullBack_to P F α' hα')
  exact ⟨InducedCategory.isoMk δ₀, InducedCategory.isoMk ε₀,
    InducedCategory.hom_ext h1, InducedCategory.hom_ext h2, InducedCategory.hom_ext h3⟩

include F in
/-- **(iii)(c)** 順方向の移送。★`𝒪^▷` の元は `.hom` でそのまま対応する。 -/
theorem istr_otriFwd {A B : Istr P} (φ : A ⟶ B) (hst : IsPreStep (istrPre P F) φ)
    (α : End A) (hα : α ∈ OTri (istrPre P F) A) :
    ∃! β : End B, β ∈ OTri (istrPre P F) B ∧ (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  obtain ⟨β₀, ⟨hβ₀m, hβ₀e⟩, hβ₀u⟩ := F.otriFwd φ.hom
    (isCoAngular_of_isotropic_dom P F A.property _) hst α.hom hα
  refine ⟨InducedCategory.homMk β₀, ⟨hβ₀m, InducedCategory.hom_ext hβ₀e⟩, ?_⟩
  rintro β ⟨hβm, hβe⟩
  exact InducedCategory.hom_ext
    (hβ₀u β.hom ⟨hβm, congrArg InducedCategory.Hom.hom hβe⟩)

include F in
/-- **(iii)(c)** 逆方向の移送。 -/
theorem istr_otriBwd {A B : Istr P} (φ : A ⟶ B) (hst : IsPreStep (istrPre P F) φ)
    (β : End B) (hβ : β ∈ OTri (istrPre P F) B) :
    ∃! α : End A, α ∈ OTri (istrPre P F) A ∧ (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  obtain ⟨α₀, ⟨hα₀m, hα₀e⟩, hα₀u⟩ := F.otriBwd φ.hom
    (isCoAngular_of_isotropic_dom P F A.property _) hst β.hom hβ
  refine ⟨InducedCategory.homMk α₀, ⟨hα₀m, InducedCategory.hom_ext hα₀e⟩, ?_⟩
  rintro α ⟨hαm, hαe⟩
  exact InducedCategory.hom_ext
    (hα₀u α.hom ⟨hαm, congrArg InducedCategory.Hom.hom hαe⟩)

include F in
/-- **(iii)(c)** `Base` にしか依らないことの移送。 -/
theorem istr_otriBase {A B : Istr P} (φ φ' : A ⟶ B)
    (hst : IsPreStep (istrPre P F) φ) (hst' : IsPreStep (istrPre P F) φ')
    (hbase : (istrPre P F).Base φ = (istrPre P F).Base φ')
    (α : End A) (hα : α ∈ OTri (istrPre P F) A)
    (β : End B) (hβ : β ∈ OTri (istrPre P F) B)
    (h : (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ) :
    (φ' ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ' :=
  InducedCategory.hom_ext
    (F.otriBase φ.hom φ'.hom (isCoAngular_of_isotropic_dom P F A.property _) hst
      (isCoAngular_of_isotropic_dom P F A.property _) hst' hbase α.hom hα β.hom hβ
      (congrArg InducedCategory.Hom.hom h))

include F in
/-- ★**base-identity 自己射を保つ**。四角形の `Base` 成分で `Base (hullMap)` を消す。 -/
theorem isotropification_baseIdentity {A : C} (e : A ⟶ A) (h : IsBaseIdentity P e) :
    IsBaseIdentity P (istrMap P F e) := by
  haveI hb : IsIso (P.Base (hullMap P F A)) := (hullMap_spec P F A).2.1.2
  have hsq := congrArg P.Base (isotropification_square P F e)
  rw [P.Base_comp, P.Base_comp] at hsq
  have he : P.Base e = 𝟙 _ := h.trans (P.Base_id A)
  show P.Base (istrMap P F e) = P.Base (𝟙 _)
  rw [P.Base_id]
  refine (cancel_epi (P.Base (hullMap P F A))).mp ?_
  rw [hsq, he, Category.id_comp, Category.comp_id]

include F in
/-- ★**(v)** の保存リスト —— `Div-identity` 自己射。

★`Φ(Base e) = id` なので、同型で共役を取っても `id` のまま。 -/
theorem isotropification_divIdentity {A : C} (e : A ⟶ A) (h : IsDivIdentity P e) :
    IsDivIdentity P (istrMap P F e) := by
  haveI hb : IsIso (P.Base (hullMap P F A)) := (hullMap_spec P F A).2.1.2
  have hsq := congrArg P.Base (isotropification_square P F e)
  rw [P.Base_comp, P.Base_comp] at hsq
  have he : Φ.map (P.Base e) = Φ.map (𝟙 (P.toElem.obj A).base) := by
    rw [← P.Base_id A]; exact h
  have hfac : P.Base (istrMap P F e)
      = inv (P.Base (hullMap P F A)) ≫ P.Base e ≫ P.Base (hullMap P F A) := by
    rw [← hsq, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  show Φ.map (P.Base (istrMap P F e)) = Φ.map (P.Base (𝟙 _))
  rw [P.Base_id, hfac]
  ext x
  show Φ.map (inv (P.Base (hullMap P F A)) ≫ P.Base e ≫ P.Base (hullMap P F A)) x
    = Φ.map (𝟙 _) x
  have he' : ∀ y, Φ.map (P.Base e) y = y := fun y => by
    rw [DFunLike.congr_fun he y, Φ.map_id]
  rw [Φ.map_comp (P.Base e ≫ P.Base (hullMap P F A)) (inv (P.Base (hullMap P F A))),
    Φ.map_comp (P.Base (hullMap P F A)) (P.Base e), he',
    ← Φ.map_comp (P.Base (hullMap P F A)) (inv (P.Base (hullMap P F A))), IsIso.inv_hom_id]

include F in
/-- ★**(v)** の保存リスト —— `base-FSM-morphism`。

★`Base` が同型で挟まれるだけなので、`fiberwise-surjective` も `mono` も保たれる。 -/
theorem isotropification_baseFSM {A B : C} (f : A ⟶ B) (h : IsBaseFSM P f) :
    IsBaseFSM P (istrMap P F f) := by
  haveI hbA : IsIso (P.Base (hullMap P F A)) := (hullMap_spec P F A).2.1.2
  haveI hbB : IsIso (P.Base (hullMap P F B)) := (hullMap_spec P F B).2.1.2
  have hsq := congrArg P.Base (isotropification_square P F f)
  rw [P.Base_comp, P.Base_comp] at hsq
  have hfac : P.Base (istrMap P F f)
      = inv (P.Base (hullMap P F A)) ≫ P.Base f ≫ P.Base (hullMap P F B) := by
    rw [← hsq, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  obtain ⟨hfs, hmono⟩ := h
  constructor
  · intro Z γ
    obtain ⟨Dd, δB, δZ, hδ⟩ := hfs (γ ≫ inv (P.Base (hullMap P F B)))
    refine ⟨Dd, δB ≫ P.Base (hullMap P F A), δZ, ?_⟩
    rw [hfac, Category.assoc, ← Category.assoc (P.Base (hullMap P F A)), IsIso.hom_inv_id,
      Category.id_comp, ← Category.assoc, hδ, Category.assoc, Category.assoc,
      IsIso.inv_hom_id, Category.comp_id]
  · rw [hfac]
    haveI := hmono
    infer_instance

include F in
/-- ★**(v)** の保存リスト —— `LB-invertible`。

★`𝒞^istr` では co-angular が自動なので、実質は isometric の移送だけ。 -/
theorem isotropification_lbInvertible {A B : C} (f : A ⟶ B) (h : IsLBInvertible P f) :
    IsLBInvertible P (istrMap P F f) :=
  ⟨isCoAngular_of_isotropic_dom P F (hullMap_spec P F A).2.2.1 _,
   (isotropification_isometric_iff P F f).mpr h.2⟩

include F in
/-- ★**(v)** の保存リスト —— `co-angular`。★**監査で「列挙どおりの形が無い」と指摘された 1 件。**

原文 (FrdI p.32):
> morphisms, base-FSM-morphisms, base-identity en-

★★**仮定は要らない** —— `𝒞^istr` の対象は isotropic なので、
そこから出る射はすべて co-angular(`isCoAngular_of_isotropic_dom`)。
★原文の「preserves co-angular morphisms」より**強い**。 -/
theorem isotropification_coAngular {A B : C} (f : A ⟶ B) :
    IsCoAngular P (istrMap P F f) :=
  isCoAngular_of_isotropic_dom P F (hullMap_spec P F A).2.2.1 _

/-! ### ★★(v) の「through which the functor `C →FΦ` factors」

原文 (FrdI p.32):
> Bistr forms a left adjoint to the inclusion functor Cistr →C, through which

★★**原文の代名詞「which」の先行詞が曖昧である**(測定として記録する)。
直前の名詞句は「the inclusion functor `Cistr →C`」だが、
★**包含関手は `𝒞` へ**入る**ので、`𝒞 → 𝔽_Φ` がそれを経由することはできない**(向きが合わない)。
★**型が通る唯一の読みは「isotropification 関手 `𝒞 → 𝒞^istr` を経由する」**である。
★我々はこちらを採る。★**原文の曖昧さであり、我々の選択である**ことを明記しておく。

★**「経由する」の中身**: `A` とその isotropic hull `A^istr` は `𝒞` では別の対象だが、
★★**`𝔽_Φ` に落とすと同型になる。** isotropic hull は
- isometric ⟹ `Div = 0`
- pre-step ⟹ `IsLinear`(`degFr = 1`)かつ base-isomorphism
なので、`𝔽_Φ` の同型判定 `isIso_iff`(base 同型 ∧ div 可逆 ∧ deg = 1)を満たす。
-/

include F in
/-- ★★**isotropic hull は `𝔽_Φ` では同型になる**。

★これが「factors」の中身である。 -/
theorem toElem_map_hullMap_isIso (A : C) : IsIso (P.toElem.map (hullMap P F A)) := by
  obtain ⟨hisom, hstep, -, -⟩ := hullMap_spec P F A
  refine (ElemFrobCat.isIso_iff _).mpr ⟨hstep.2, ?_, hstep.1⟩
  show IsAddUnit (P.Div (hullMap P F A))
  rw [show P.Div (hullMap P F A) = 0 from hisom]
  exact isAddUnit_zero

include F in
/-- ★★**`Proposition 1.9, (v)` の「through which the functor `C →FΦ` factors」**。

原文 (FrdI p.32):
> the functor C →FΦ factors.

★**`𝒞 → 𝔽_Φ` は isotropification を経由する**(自然同型を除いて)。
★成分は `toElem_map_hullMap_isIso`、自然性は `isotropification_square` そのもの。 -/
noncomputable def isotropificationFactorIso :
    P.toElem ≅ isotropification P F ⋙ (istrPre P F).toElem :=
  NatIso.ofComponents
    (fun A => @asIso _ _ _ _ (P.toElem.map (hullMap P F A)) (toElem_map_hullMap_isIso P F A))
    (fun {A B} f => by
      show P.toElem.map f ≫ P.toElem.map (hullMap P F B)
        = P.toElem.map (hullMap P F A) ≫ P.toElem.map (istrMap P F f)
      rw [← P.toElem.map_comp, ← P.toElem.map_comp, isotropification_square P F f])

include F in
/-- ★**(v)** —— **isotropification 関手の `𝒞^istr` への制限は恒等関手と同型**。

原文 (FrdI p.32):
> morphic to the identity functor. Finally, the isotropification functor preserves

★成分は `hullMap_isIso`(isotropic な対象では hull への射が同型)、
自然性は `isotropification_square` そのもの。 -/
noncomputable def isotropificationRestrictIso :
    (isotropicProp P).ι ⋙ isotropification P F ≅ 𝟭 (Istr P) :=
  NatIso.ofComponents
    (fun X => (InducedCategory.isoMk
      (@asIso _ _ _ _ (hullMap P F X.obj) (hullMap_isIso P F X.obj X.property))).symm)
    (fun {X Y} g => by
      haveI hX : IsIso (hullMap P F X.obj) := hullMap_isIso P F X.obj X.property
      haveI hY : IsIso (hullMap P F Y.obj) := hullMap_isIso P F Y.obj Y.property
      refine InducedCategory.hom_ext ?_
      show istrMap P F g.hom ≫ inv (hullMap P F Y.obj)
        = inv (hullMap P F X.obj) ≫ g.hom
      rw [IsIso.comp_inv_eq, Category.assoc, IsIso.eq_inv_comp]
      exact isotropification_square P F g.hom)

include F in
/-- **(i)(a)** の移送(★**後ろ向き** —— isotropification が要る)。

`𝒞` で得た Frobenius-trivial な対象 `A₀` の isotropic hull を取る。
* Frobenius-trivial の `ζ` は **isotropification が関手であること**でそのまま運べる
* 底の同型は `Base (hullMap)` が同型であることで繋ぐ -/
theorem istr_baseSurj (Y : D) :
    ∃ A : Istr P, IsFrobeniusTrivial (istrPre P F) A ∧
      Nonempty (((istrPre P F).toElem.obj A).base ≅ Y) := by
  obtain ⟨A₀, ⟨ζ, hdeg, hprop⟩, ⟨e⟩⟩ := F.baseSurj Y
  haveI hb : IsIso (P.Base (hullMap P F A₀)) := (hullMap_spec P F A₀).2.1.2
  refine ⟨hullIstr P F A₀, ⟨⟨⟨fun n => (isotropification P F).map (ζ n), ?_⟩, ?_⟩, ?_, ?_⟩,
    ⟨(asIso (P.Base (hullMap P F A₀))).symm ≪≫ e⟩⟩
  · show (isotropification P F).map (ζ 1) = 𝟙 _
    rw [show ζ 1 = 𝟙 A₀ from ζ.map_one]
    exact (isotropification P F).map_id _
  · intro m n
    show (isotropification P F).map (ζ (m * n))
      = (isotropification P F).map (ζ n) ≫ (isotropification P F).map (ζ m)
    rw [ζ.map_mul]
    exact (isotropification P F).map_comp (ζ n) (ζ m)
  · intro n
    show P.degFr (istrMap P F (ζ n)) = n
    rw [isotropification_degFr]
    exact hdeg n
  · intro n
    exact ⟨isotropification_baseIdentity P F (ζ n) (hprop n).1,
      (istr_frobType_iff P F _).mpr (isotropification_frobType P F (ζ n) (hprop n).2)⟩

include F in
/-- **(i)(b)** の移送(★**後ろ向き** —— isotropification が要る)。

`𝒞` で得た span `X₀ ⟶ A.obj`, `X₀ ⟶ B.obj` を isotropification で運び、
`A`, `B` が既に isotropic なので `hullMap` の逆で戻す。

★★**手順2(`inv` を書かない)の例外条件**: **主張そのものが
`@inv _ _ _ _ (Base φ) hφ.2` を含む**ので避けられない。
`IsIso.eq_inv_comp` を**インスタンスを明示して**当てる形で扱う。
★一方、証明の中で使う逆射はすべて `IsIso.out` から明示的に取り、`inv` は書かない。
★**手順は固定するものではなく、反例が出たら条件を精密化するもの** という形の一例。 -/
theorem istr_preStepSpan (A B : Istr P)
    (α : ((istrPre P F).toElem.obj A).base ⟶ ((istrPre P F).toElem.obj B).base) (hα : IsIso α) :
    ∃ (X : Istr P) (φ : X ⟶ A) (ψ : X ⟶ B) (hφ : IsPreStep (istrPre P F) φ),
      IsPreStep (istrPre P F) ψ ∧
        α = @inv _ _ _ _ ((istrPre P F).Base φ) hφ.2 ≫ (istrPre P F).Base ψ := by
  obtain ⟨X₀, φ₀, ψ₀, hφ₀, hψ₀, heq⟩ := F.preStepSpan A.obj B.obj α hα
  haveI hmA : IsIso (hullMap P F A.obj) := hullMap_isIso P F A.obj A.property
  haveI hmB : IsIso (hullMap P F B.obj) := hullMap_isIso P F B.obj B.property
  obtain ⟨wA, hwA1, -⟩ := hmA.out
  obtain ⟨wB, hwB1, -⟩ := hmB.out
  haveI hiA : IsIso wA := ⟨hullMap P F A.obj, by
      refine (cancel_epi (hullMap P F A.obj)).mp ?_
      rw [← Category.assoc, hwA1, Category.id_comp, Category.comp_id], hwA1⟩
  haveI hiB : IsIso wB := ⟨hullMap P F B.obj, by
      refine (cancel_epi (hullMap P F B.obj)).mp ?_
      rw [← Category.assoc, hwB1, Category.id_comp, Category.comp_id], hwB1⟩
  haveI hXb : IsIso (P.Base (hullMap P F X₀)) := (hullMap_spec P F X₀).2.1.2
  obtain ⟨v, hv1, hv2⟩ := hXb.out
  -- ★`Base` の計算: hull を通しても `v` を挟むだけ
  have key : ∀ (Z : C) (f : X₀ ⟶ Z) (w : hullObj P F Z ⟶ Z),
      hullMap P F Z ≫ w = 𝟙 Z → P.Base (istrMap P F f ≫ w) = v ≫ P.Base f := by
    intro Z f w hw
    refine (cancel_epi (P.Base (hullMap P F X₀))).mp ?_
    have hsq := congrArg P.Base (isotropification_square P F f)
    rw [P.Base_comp, P.Base_comp] at hsq
    calc P.Base (hullMap P F X₀) ≫ P.Base (istrMap P F f ≫ w)
        = P.Base (hullMap P F X₀) ≫ (P.Base (istrMap P F f) ≫ P.Base w) := by
          rw [P.Base_comp]
      _ = (P.Base (hullMap P F X₀) ≫ P.Base (istrMap P F f)) ≫ P.Base w :=
          (Category.assoc _ _ _).symm
      _ = (P.Base f ≫ P.Base (hullMap P F Z)) ≫ P.Base w := by rw [hsq]
      _ = P.Base f ≫ (P.Base (hullMap P F Z) ≫ P.Base w) := Category.assoc _ _ _
      _ = P.Base f ≫ P.Base (hullMap P F Z ≫ w) := by rw [P.Base_comp]
      _ = P.Base f ≫ P.Base (𝟙 Z) := by rw [hw]
      _ = P.Base f := by rw [P.Base_id, Category.comp_id]
      _ = (P.Base (hullMap P F X₀) ≫ v) ≫ P.Base f := by rw [hv1, Category.id_comp]
      _ = P.Base (hullMap P F X₀) ≫ (v ≫ P.Base f) := Category.assoc _ _ _
  have hpsA : IsPreStep P (istrMap P F φ₀ ≫ wA) :=
    IsPreStep.comp P ((isotropification_preStep_iff P F φ₀).mpr hφ₀)
      (isPreStep_of_isIso P wA)
  have hpsB : IsPreStep P (istrMap P F ψ₀ ≫ wB) :=
    IsPreStep.comp P ((isotropification_preStep_iff P F ψ₀).mpr hψ₀)
      (isPreStep_of_isIso P wB)
  refine ⟨hullIstr P F X₀, InducedCategory.homMk (istrMap P F φ₀ ≫ wA),
    InducedCategory.homMk (istrMap P F ψ₀ ≫ wB), hpsA, hpsB, ?_⟩
  -- `heq` を `Base φ₀ ≫ α = Base ψ₀` に直す
  have heq' : P.Base φ₀ ≫ α = P.Base ψ₀ :=
    (@IsIso.eq_inv_comp _ _ _ _ _ (P.Base φ₀) hφ₀.2 _ _).mp heq
  refine (@IsIso.eq_inv_comp _ _ _ _ _ ((istrPre P F).Base
    (InducedCategory.homMk (istrMap P F φ₀ ≫ wA) : hullIstr P F X₀ ⟶ A)) hpsA.2 _ _).mpr ?_
  show P.Base (istrMap P F φ₀ ≫ wA) ≫ α = P.Base (istrMap P F ψ₀ ≫ wB)
  rw [key A.obj φ₀ wA hwA1, key B.obj ψ₀ wB hwB1, Category.assoc, heq']
  rfl

/-! ### ★**(i)(c) の圏同値の移送** —— 21 条の最後 -/

include F in
/-- 補助: `(𝒞^istr)^pl-bk` の `A` 上の対象を `𝒞^pl-bk` の `A.obj` 上へ運ぶ。

★ここで **`istr_isPullBack_to`(難しい向き)** を使う。 -/
def istrPlBkToC (A : Istr P) (Z : Over (⟨A⟩ : PlBk (istrPre P F))) :
    Over (⟨A.obj⟩ : PlBk P) :=
  Over.mk (⟨Z.hom.hom.hom, istr_isPullBack_to P F Z.hom.hom Z.hom.property⟩ :
    (⟨Z.left.obj.obj⟩ : PlBk P) ⟶ (⟨A.obj⟩ : PlBk P))

include F in
/-- 補助: 射のほうの運搬。 -/
def istrPlBkToCMap {A : Istr P} {Z W : Over (⟨A⟩ : PlBk (istrPre P F))} (h : Z ⟶ W) :
    istrPlBkToC P F A Z ⟶ istrPlBkToC P F A W :=
  Over.homMk (⟨h.left.hom.hom, istr_isPullBack_to P F h.left.hom h.left.property⟩ :
      (⟨Z.left.obj.obj⟩ : PlBk P) ⟶ (⟨W.left.obj.obj⟩ : PlBk P))
    (by
      have hw : h.left.hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w h)
      have hw2 : h.left.hom.hom ≫ W.hom.hom.hom = Z.hom.hom.hom :=
        congrArg InducedCategory.Hom.hom hw
      exact InducedWideCategory.Hom.ext hw2)

include F in
/-- **(i)(c)** の移送。★`𝒞^istr` は `𝒞` の**充満**部分圏なので、
`(𝒞^istr)^pl-bk` の `A` 上のスライスは `𝒞^pl-bk` の `A.obj` 上のスライスと
**圏として一致する** —— 対象が一致するのは
`isotropic_dom_of_pullBack`(isotropic への pull-back の始域は isotropic)による。

★充満性・忠実性・本質的全射性の**3つとも** `F.plBkEquiv A.obj` から来る。 -/
theorem istr_plBkEquiv (A : Istr P) :
    (plBkOverFunctor (istrPre P F) A).IsEquivalence := by
  haveI := F.plBkEquiv A.obj
  haveI hfaith : (plBkOverFunctor (istrPre P F) A).Faithful := by
    constructor
    intro Z W f g hfg
    have hb : (istrPre P F).Base f.left.hom = (istrPre P F).Base g.left.hom :=
      congrArg CommaMorphism.left hfg
    have hmap : (plBkOverFunctor P A.obj).map (istrPlBkToCMap P F f)
        = (plBkOverFunctor P A.obj).map (istrPlBkToCMap P F g) :=
      Over.OverMorphism.ext hb
    have h2 := (plBkOverFunctor P A.obj).map_injective hmap
    have h3 : f.left.hom.hom = g.left.hom.hom :=
      congrArg (fun x => InducedWideCategory.Hom.hom (CommaMorphism.left x)) h2
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext (InducedCategory.hom_ext h3))
  haveI hfull : (plBkOverFunctor (istrPre P F) A).Full := by
    constructor
    intro Z W h
    obtain ⟨f', hf'⟩ := (plBkOverFunctor P A.obj).map_surjective
      (show (plBkOverFunctor P A.obj).obj (istrPlBkToC P F A Z) ⟶
          (plBkOverFunctor P A.obj).obj (istrPlBkToC P F A W) from h)
    have hw : f'.left.hom ≫ W.hom.hom.hom = Z.hom.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f')
    refine ⟨Over.homMk (⟨InducedCategory.homMk f'.left.hom,
        istr_isPullBack_of P F (InducedCategory.homMk f'.left.hom) f'.left.property⟩ :
        Z.left ⟶ W.left)
      (InducedWideCategory.Hom.ext (InducedCategory.hom_ext hw)), ?_⟩
    exact Over.OverMorphism.ext (congrArg CommaMorphism.left hf')
  haveI hess : (plBkOverFunctor (istrPre P F) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨Z', hZ'⟩ : ∃ Z' : Over (⟨A.obj⟩ : PlBk P),
        Z' = (plBkOverFunctor P A.obj).objPreimage Y := ⟨_, rfl⟩
    have hiZ : (plBkOverFunctor P A.obj).obj Z' ≅ Y := by
      rw [hZ']; exact (plBkOverFunctor P A.obj).objObjPreimageIso Y
    have hiso : IsIsotropic P Z'.left.obj :=
      isotropic_dom_of_pullBack P F Z'.hom.hom Z'.hom.property A.property
    refine ⟨Over.mk (⟨(InducedCategory.homMk Z'.hom.hom :
          (⟨Z'.left.obj, hiso⟩ : Istr P) ⟶ A),
        istr_isPullBack_of P F _ Z'.hom.property⟩ :
      (⟨(⟨Z'.left.obj, hiso⟩ : Istr P)⟩ : PlBk (istrPre P F)) ⟶ (⟨A⟩ : PlBk (istrPre P F))), ?_⟩
    exact ⟨hiZ⟩
  exact ⟨hfaith, hfull, hess⟩

include F in
/-- ★★★**`𝒞^istr` は `Definition 1.3` の core 21 条をすべて満たす**。

原文 (FrdI p.32):
> (v) Cistr [equipped with the restriction to C of the given functor C →FΦ] is a

★**21 条のうち 19 条は「前向き」で自動**((vii)(b) が `𝒞^istr` を閉じるから)、
**2 条だけが「後ろ向き」**で isotropification を要した(`baseSurj` / `preStepSpan`)。
`plBkEquiv` は両向きが要った(`istr_isPullBack_to` が随伴を使う)。 -/
theorem istr_frobenioidCore : FrobenioidCore (istrPre P F) where
  baseSurj := istr_baseSurj P F
  preStepSpan := istr_preStepSpan P F
  plBkEquiv := istr_plBkEquiv P F
  frobDegSurj := istr_frobDegSurj P F
  frobDegUniq := istr_frobDegUniq P F
  coAngularComp := istr_coAngularComp P F
  coAngularOfPreStep := fun α hca hst φ => istr_coAngularOfPreStep P F α hca hst φ
  otriFwd := fun φ _ hst α hα => istr_otriFwd P F φ hst α hα
  otriBwd := fun φ _ hst β hβ => istr_otriBwd P F φ hst β hβ
  otriBase := fun φ φ' _ hst _ hst' hbase α hα β hβ h =>
    istr_otriBase P F φ φ' hst hst' hbase α hα β hβ h
  arbFactor := istr_arbFactor P F
  arbFactorUniq := istr_arbFactorUniq P F
  pullBackLB := fun α h => istr_pullBackLB P F α h
  preStepMono := istr_preStepMono P F
  preStepFactor := istr_preStepFactor P F
  preStepFactorUniq := istr_preStepFactorUniq P F
  preStepFactor' := istr_preStepFactor' P F
  preStepFactorUniq' := istr_preStepFactorUniq' P F
  faithfulUpToUnits := istr_faithfulUpToUnits P F
  isotropicHullExists := istr_isotropicHullExists P F
  isotropicClosed := istr_isotropicClosed P F

/-! ### ★**(iii)(d) の圏同値2本の移送**

★`plBkEquiv` と**同じ形**である: 忠実性だけは `𝒞^istr` の中で直接出て
(全射性/単射性 —— `𝒞^istr` も totally epimorphic、pre-step も mono)、
**充満性と本質的全射性は `𝒞` 側の同値から引く**。
違うのは「対象が `𝒞^istr` に残ること」を言う補題だけ:

* コスライス側は **`isotropicClosed`**(前向き、自明)、
* スライス側は **`Proposition 1.9, (iv)`**(co-angular pre-step は co-angular linear)。
-/

section CoaPre

variable [(coaPreProp P).IsMultiplicative] [(coaPreProp (istrPre P F)).IsMultiplicative]

include F in
/-- 補助: `_A(𝒞^istr,coa-pre)` の対象を `_{A.obj}(𝒞^coa-pre)` へ。 -/
def istrCoaPreUnder (A : Istr P)
    (Z : Under (⟨A⟩ : WideSubcategory (coaPreProp (istrPre P F)))) :
    Under (⟨A.obj⟩ : WideSubcategory (coaPreProp P)) :=
  Under.mk (⟨Z.hom.hom.hom,
      isCoAngular_of_isotropic_dom P F A.property _, Z.hom.property.2⟩ :
    (⟨A.obj⟩ : WideSubcategory (coaPreProp P)) ⟶
      (⟨Z.right.obj.obj⟩ : WideSubcategory (coaPreProp P)))

include F in
/-- 補助: `(𝒞^istr,coa-pre)_A` の対象を `(𝒞^coa-pre)_{A.obj}` へ。 -/
def istrCoaPreOver (A : Istr P)
    (Z : Over (⟨A⟩ : WideSubcategory (coaPreProp (istrPre P F)))) :
    Over (⟨A.obj⟩ : WideSubcategory (coaPreProp P)) :=
  Over.mk (⟨Z.hom.hom.hom,
      isCoAngular_of_isotropic_dom P F Z.left.obj.property _, Z.hom.property.2⟩ :
    (⟨Z.left.obj.obj⟩ : WideSubcategory (coaPreProp P)) ⟶
      (⟨A.obj⟩ : WideSubcategory (coaPreProp P)))

include F in
/-- **(iii)(d)** コスライス側の移送。 -/
theorem istr_coaPreUnderEquiv
    (hC : ∀ A : C, (coaPreUnderFunctor P A).IsEquivalence) (A : Istr P) :
    (coaPreUnderFunctor (istrPre P F) A).IsEquivalence := by
  haveI := hC A.obj
  haveI hfaith : (coaPreUnderFunctor (istrPre P F) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : Z.hom.hom ≫ f.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f)
    have h2 : Z.hom.hom ≫ g.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w g)
    haveI : Epi Z.hom.hom := (istrPre P F).totEpiC _ _ _
    exact Under.UnderMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_epi Z.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreUnderFunctor (istrPre P F) A).Full := by
    constructor
    intro Z W h
    obtain ⟨f', -⟩ := (coaPreUnderFunctor P A.obj).map_surjective
      (show (coaPreUnderFunctor P A.obj).obj (istrCoaPreUnder P F A Z) ⟶
          (coaPreUnderFunctor P A.obj).obj (istrCoaPreUnder P F A W) from h)
    have hw : Z.hom.hom.hom ≫ f'.right.hom = W.hom.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f')
    refine ⟨Under.homMk (⟨(InducedCategory.homMk f'.right.hom :
          Z.right.obj ⟶ W.right.obj),
        ⟨istr_coAngular P F _, f'.right.property.2⟩⟩ : Z.right ⟶ W.right)
      (InducedWideCategory.Hom.ext (InducedCategory.hom_ext hw)), Subsingleton.elim _ _⟩
  haveI hess : (coaPreUnderFunctor (istrPre P F) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨Z', hZ'⟩ : ∃ Z' : Under (⟨A.obj⟩ : WideSubcategory (coaPreProp P)),
        Z' = (coaPreUnderFunctor P A.obj).objPreimage Y := ⟨_, rfl⟩
    have hiZ : (coaPreUnderFunctor P A.obj).obj Z' ≅ Y := by
      rw [hZ']; exact (coaPreUnderFunctor P A.obj).objObjPreimageIso Y
    have hiso : IsIsotropic P Z'.right.obj := F.isotropicClosed Z'.hom.hom A.property
    exact ⟨Under.mk (⟨(InducedCategory.homMk Z'.hom.hom :
          A ⟶ (⟨Z'.right.obj, hiso⟩ : Istr P)),
        ⟨istr_coAngular P F _, Z'.hom.property.2⟩⟩ :
      (⟨A⟩ : WideSubcategory (coaPreProp (istrPre P F))) ⟶
        (⟨(⟨Z'.right.obj, hiso⟩ : Istr P)⟩ :
          WideSubcategory (coaPreProp (istrPre P F)))), ⟨hiZ⟩⟩
  exact ⟨hfaith, hfull, hess⟩

include F in
/-- **(iii)(d)** スライス側の移送。★対象が `𝒞^istr` に残ることに
**`Proposition 1.9, (iv)` を使う** —— co-angular pre-step は co-angular linear なので、
終域が isotropic なら始域も isotropic。 -/
theorem istr_coaPreOverEquiv
    (hC : ∀ A : C, (coaPreOverFunctor P A).IsEquivalence) (A : Istr P) :
    (coaPreOverFunctor (istrPre P F) A).IsEquivalence := by
  haveI := hC A.obj
  haveI hfaith : (coaPreOverFunctor (istrPre P F) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have h2 : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    haveI : Mono W.hom.hom := istr_preStepMono P F _ W.hom.property.2
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_mono W.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreOverFunctor (istrPre P F) A).Full := by
    constructor
    intro Z W h
    obtain ⟨f', -⟩ := (coaPreOverFunctor P A.obj).map_surjective
      (show (coaPreOverFunctor P A.obj).obj (istrCoaPreOver P F A Z) ⟶
          (coaPreOverFunctor P A.obj).obj (istrCoaPreOver P F A W) from h)
    have hw : f'.left.hom ≫ W.hom.hom.hom = Z.hom.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f')
    refine ⟨Over.homMk (⟨(InducedCategory.homMk f'.left.hom :
          Z.left.obj ⟶ W.left.obj),
        ⟨istr_coAngular P F _, f'.left.property.2⟩⟩ : Z.left ⟶ W.left)
      (InducedWideCategory.Hom.ext (InducedCategory.hom_ext hw)), Subsingleton.elim _ _⟩
  haveI hess : (coaPreOverFunctor (istrPre P F) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨Z', hZ'⟩ : ∃ Z' : Over (⟨A.obj⟩ : WideSubcategory (coaPreProp P)),
        Z' = (coaPreOverFunctor P A.obj).objPreimage Y := ⟨_, rfl⟩
    have hiZ : (coaPreOverFunctor P A.obj).obj Z' ≅ Y := by
      rw [hZ']; exact (coaPreOverFunctor P A.obj).objObjPreimageIso Y
    have hiso : IsIsotropic P Z'.left.obj :=
      (prop_1_9_iv P F Z'.hom.hom Z'.hom.property.1 Z'.hom.property.2.1).mpr A.property
    exact ⟨Over.mk (⟨(InducedCategory.homMk Z'.hom.hom :
          (⟨Z'.left.obj, hiso⟩ : Istr P) ⟶ A),
        ⟨istr_coAngular P F _, Z'.hom.property.2⟩⟩ :
      (⟨(⟨Z'.left.obj, hiso⟩ : Istr P)⟩ : WideSubcategory (coaPreProp (istrPre P F))) ⟶
        (⟨A⟩ : WideSubcategory (coaPreProp (istrPre P F)))), ⟨hiZ⟩⟩
  exact ⟨hfaith, hfull, hess⟩

end CoaPre

/-- ★★★**`Proposition 1.9, (v)`** —— `𝒞^istr` は Frobenioid である。 -/
theorem istr_frobenioid (G : Frobenioid P) : Frobenioid (istrPre P F) := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := coaPreProp_isMultiplicative (istrPre P F)
    (istr_frobenioidCore P F).coAngularComp
  exact ⟨istr_frobenioidCore P F,
    istr_coaPreUnderEquiv P F G.coaPreUnderEquiv,
    istr_coaPreOverEquiv P F G.coaPreOverEquiv⟩

end Istr

/-! ## ★第4段 —— (ii)(iii) の梱包: `φ_*` と `φ^*` を `Functor` にする

原文 (FrdI p.31):
> (ii) Any base-isomorphism φ : A →B of C induces a functor [well-defined

原文 (FrdI p.31):
> (iii) Any pull-back morphism φ : A →B of C induces a functor [well-

★★**測定**: 関手の**対象の割り当てにだけ選択が要り**、
**射の割り当ては一意**(`imtrPre_hom_uniq`)、
**関手則(`map_id` / `map_comp`)は無料**である ——
`𝒞^imtr-pre_A` の hom 集合が**高々1元**だからである。
★原文の「[well-defined up to isomorphism]」という但し書きは
**対象の割り当てだけに掛かる**。
-/

section ImtrPreFunctor

variable (F : FrobenioidCore P)

include F in
/-- ★★**`𝒞^imtr-pre_A` は thin**(hom 集合が高々1元)。

★`Definition 1.3, (v), (a)`(pre-step は mono)だけから出る。
これにより `φ_*` / `φ^*` の**忠実性と関手則が無料**になる。

★**手順3**(移送・補助補題を書くときは最初に `include F in`)—— ここでも同じ。 -/
theorem imtrPreOver_hom_subsingleton {A : C} {Z W : Over (⟨A⟩ : ImtrPre P)} (f g : Z ⟶ W) :
    f = g :=
  Over.OverMorphism.ext (InducedWideCategory.Hom.ext
    (imtrPre_hom_uniq P F Z.hom.hom W.hom.hom W.hom.property.2 _ _
      (congrArg InducedWideCategory.Hom.hom (Over.w f))
      (congrArg InducedWideCategory.Hom.hom (Over.w g))))

/-! ### `φ_*`(base-isomorphism から) -/

/-- ★分解の左因子(co-angular base-isomorphism)。 -/
noncomputable def pushMid {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Cc : C} (ε : Cc ⟶ A) (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    Cc ⟶ pushObj P F φ hφ ε hεi hεs :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose

theorem pushFac {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Cc : C} (ε : Cc ⟶ A) (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    ε ≫ φ = pushMid P F φ hφ ε hεi hεs ≫ pushHom P F φ hφ ε hεi hεs :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose_spec.choose_spec.1

theorem pushMid_spec {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Cc : C} (ε : Cc ⟶ A) (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    IsCoAngular P (pushMid P F φ hφ ε hεi hεs) ∧
      IsBaseIsomorphism P (pushMid P F φ hφ ε hεi hεs) :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose_spec.choose_spec.2.1

/-- ★`φ_*` の射への割り当て。★**選択は要らない**が、`∃` から取るので
`.choose` を使う —— 一意性は `imtrPre_hom_uniq` にある。 -/
noncomputable def pushMap {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Z W : Over (⟨A⟩ : ImtrPre P)} (f : Z ⟶ W) :
    pushObj P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2 ⟶
      pushObj P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2 :=
  (prop_1_9_ii_hom P F φ
    Z.hom.hom (pushMid P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
      (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    W.hom.hom (pushMid P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
      (pushHom P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
    (pushFac P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    (pushMid_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushMid_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushFac P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
    (pushMid_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).2
    (pushHom_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).1
    (pushHom_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).2
    f.left.hom f.left.property.2
    (congrArg InducedWideCategory.Hom.hom (Over.w f))).choose

theorem pushMap_spec {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Z W : Over (⟨A⟩ : ImtrPre P)} (f : Z ⟶ W) :
    (IsIsometric P (pushMap P F φ hφ f) ∧ IsPreStep P (pushMap P F φ hφ f)) ∧
      pushMap P F φ hφ f ≫
          pushHom P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2 =
        pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2 :=
  (prop_1_9_ii_hom P F φ
    Z.hom.hom (pushMid P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
      (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    W.hom.hom (pushMid P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
      (pushHom P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
    (pushFac P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    (pushMid_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushMid_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushFac P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
    (pushMid_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).2
    (pushHom_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).1
    (pushHom_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).2
    f.left.hom f.left.property.2
    (congrArg InducedWideCategory.Hom.hom (Over.w f))).choose_spec

/-- ★★**`φ_* : 𝒞^imtr-pre_A ⥤ 𝒞^imtr-pre_B`**(`φ` は base-isomorphism)。 -/
noncomputable def pushFunctor {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ) :
    Over (⟨A⟩ : ImtrPre P) ⥤ Over (⟨B⟩ : ImtrPre P) where
  obj Z := Over.mk (⟨pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2,
      pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ :
    (⟨pushObj P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ : ImtrPre P) ⟶
      (⟨B⟩ : ImtrPre P))
  map {Z W} f := Over.homMk (⟨pushMap P F φ hφ f, (pushMap_spec P F φ hφ f).1⟩ :
      (⟨pushObj P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ : ImtrPre P) ⟶
        (⟨pushObj P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2⟩ : ImtrPre P))
    (InducedWideCategory.Hom.ext (pushMap_spec P F φ hφ f).2)
  map_id _ := imtrPreOver_hom_subsingleton P F _ _
  map_comp _ _ := imtrPreOver_hom_subsingleton P F _ _

/-! ### `φ^*`(pull-back morphism から) -/

noncomputable def pullObj {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) : C :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose

noncomputable def pullHom {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    pullObj P F φ hpb δ hδi hδs ⟶ A :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose

/-- ★四角形の下辺(pull-back)。 -/
noncomputable def pullPsi {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    pullObj P F φ hpb δ hδi hδs ⟶ Dd :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose

theorem pullFac {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    pullHom P F φ hpb δ hδi hδs ≫ φ = pullPsi P F φ hpb δ hδi hδs ≫ δ :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose_spec.1

theorem pullPsi_spec {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    IsPullBack P (pullPsi P F φ hpb δ hδi hδs) :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose_spec.2.1

theorem pullHom_spec {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    IsIsometric P (pullHom P F φ hpb δ hδi hδs) ∧
      IsPreStep P (pullHom P F φ hpb δ hδi hδs) :=
  ⟨(prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose_spec.2.2.1,
   (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose_spec.2.2.2⟩

/-- ★`φ^*` の射への割り当て。 -/
noncomputable def pullMap {A B : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    {Z W : Over (⟨B⟩ : ImtrPre P)} (f : Z ⟶ W) :
    pullObj P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ⟶
      pullObj P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2 :=
  (prop_1_9_iii_hom P φ hpb W.hom.hom W.hom.property.2
    (pullHom P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    (pullPsi P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ f.left.hom)
    (pullHom P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (pullPsi P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (by
      have hw : (Over.Hom.left f).hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w f)
      rw [pullFac P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2, Category.assoc, hw])
    (pullFac P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (pullHom_spec P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2).2
    (pullPsi_spec P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)).choose

theorem pullMap_tri {A B : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    {Z W : Over (⟨B⟩ : ImtrPre P)} (f : Z ⟶ W) :
    pullMap P F φ hpb f ≫
        pullHom P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2 =
      pullHom P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2 :=
  (prop_1_9_iii_hom P φ hpb W.hom.hom W.hom.property.2
    (pullHom P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    (pullPsi P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ f.left.hom)
    (pullHom P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (pullPsi P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (by
      have hw : (Over.Hom.left f).hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w f)
      rw [pullFac P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2, Category.assoc, hw])
    (pullFac P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (pullHom_spec P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2).2
    (pullPsi_spec P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)).choose_spec.1

/-- ★`pullMap` が本当に isometric pre-step であることは
**`Proposition 1.7, (v)`**(合成が属せば因子も属する)から出る ——
`prop_1_9_iii_hom` はそこまでは言わない。 -/
theorem pullMap_spec {A B : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    {Z W : Over (⟨B⟩ : ImtrPre P)} (f : Z ⟶ W) :
    IsIsometric P (pullMap P F φ hpb f) ∧ IsPreStep P (pullMap P F φ hpb f) := by
  have ht := pullMap_tri P F φ hpb f
  have h1 := (pullHom_spec P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
  have h2 := (pullHom_spec P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
  rw [← ht] at h1 h2
  exact ⟨(prop_1_7_v_isometric P _ _ h1).1, (prop_1_7_v_preStep P _ _ h2).1⟩

/-- ★★**`φ^* : 𝒞^imtr-pre_B ⥤ 𝒞^imtr-pre_A`**(`φ` は pull-back morphism)。 -/
noncomputable def pullFunctor {A B : C} (φ : A ⟶ B) (hpb : IsPullBack P φ) :
    Over (⟨B⟩ : ImtrPre P) ⥤ Over (⟨A⟩ : ImtrPre P) where
  obj Z := Over.mk (⟨pullHom P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2,
      pullHom_spec P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ :
    (⟨pullObj P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ : ImtrPre P) ⟶
      (⟨A⟩ : ImtrPre P))
  map {Z W} f := Over.homMk (⟨pullMap P F φ hpb f, pullMap_spec P F φ hpb f⟩ :
      (⟨pullObj P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ : ImtrPre P) ⟶
        (⟨pullObj P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2⟩ : ImtrPre P))
    (InducedWideCategory.Hom.ext (pullMap_tri P F φ hpb f))
  map_id _ := imtrPreOver_hom_subsingleton P F _ _
  map_comp _ _ := imtrPreOver_hom_subsingleton P F _ _

/-- **`𝒪^×(A)^imtr-pre ⊆ 𝒪^×(A)`** —— `u_*` が恒等関手と同型になる `u` の全体。

原文 (FrdI p.31):
> the subgroup of v ∈O×(A) for which vimtr-pre is the identity.

★`u ∈ 𝒪^×(A)` は同型なので base-isomorphism、したがって `u_*` が定義できる。 -/
def OTimesImtrPre (A : C) : Set (End A) :=
  {u | u ∈ OTimes P A ∧ ∃ hb : IsBaseIsomorphism P (u : A ⟶ A),
    Nonempty (pushFunctor P F (u : A ⟶ A) hb ≅ 𝟭 (Over (⟨A⟩ : ImtrPre P)))}

theorem otimesImtrPre_subset (A : C) : OTimesImtrPre P F A ⊆ OTimes P A :=
  fun _ h => h.1

/-! ### ★★`𝒪^×(A)^{imtr-pre}` が**部分群**であること

原文 (FrdI p.31):
> the subgroup of v ∈O×(A) for which vimtr-pre is the identity.

★★**原文は「部分群」と述べている。** 上の `Set` だけでは主張を満たさない
(2026-08-16 の監査で判明)。単位元・積・逆元の閉性を示す。

★**楽になる事実**: `Over (⟨A⟩ : ImtrPre P)` の hom は subsingleton
(`imtrPreOver_hom_subsingleton`)。★**したがって自然性は自動**で、
対象ごとの同型さえ作れば関手の同型が得られる。 -/

/-- ★**`ImtrPre P` の同型を、`𝒞` の同型から作る**。

★同型は isometric(`isIsometric_of_isIso`)かつ pre-step(`isPreStep_of_isIso`)。 -/
def imtrPreIsoOfIso {X Y : C} (δ : X ≅ Y) : (⟨X⟩ : ImtrPre P) ≅ (⟨Y⟩ : ImtrPre P) where
  hom := ⟨δ.hom, ⟨isIsometric_of_isIso P δ.hom, isPreStep_of_isIso P δ.hom⟩⟩
  inv := ⟨δ.inv, ⟨isIsometric_of_isIso P δ.inv, isPreStep_of_isIso P δ.inv⟩⟩
  hom_inv_id := InducedWideCategory.Hom.ext (by simp)
  inv_hom_id := InducedWideCategory.Hom.ext (by simp)

include F in
/-- ★★**`(𝟙 A)_* ≅ 𝟭`** —— `𝒪^×(A)^{imtr-pre}` が**単位元を含む**ことの中身。

★`ε` 自身が isometric pre-step なので、`ε = 𝟙 ≫ ε` が第2の分解を与える。
★**分解の一意性**(`prop_1_9_i_uniq`)が対象ごとの同型を出し、
自然性は hom の subsingleton 性から自動。 -/
theorem pushId_uniq (A : C) (hb : IsBaseIsomorphism P (𝟙 A))
    (Z : Over (⟨A⟩ : ImtrPre P)) :
    ∃ δ : pushObj P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≅ Z.left.obj,
      Z.hom.hom
        = δ.inv ≫ pushHom P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ∧
      (𝟙 Z.left.obj : Z.left.obj ⟶ Z.left.obj)
        = pushMid P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ δ.hom :=
  prop_1_9_i_uniq P F _ _ _ _ (𝟙 Z.left.obj) Z.hom.hom
    (by rw [← pushFac P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2]; simp)
    (pushMid_spec P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushMid_spec P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushHom_spec P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushHom_spec P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (isCoAngular_id P _) (isBaseIsomorphism_of_isIso P _)
    Z.hom.property.1 Z.hom.property.2

include F in
/-- ★★**`(𝟙 A)_* ≅ 𝟭`**。★`∃` からデータを取るので `choose` を経由する。 -/
noncomputable def pushFunctorIdIso (A : C) (hb : IsBaseIsomorphism P (𝟙 A)) :
    pushFunctor P F (𝟙 A) hb ≅ 𝟭 (Over (⟨A⟩ : ImtrPre P)) :=
  NatIso.ofComponents
    (fun Z =>
      Over.isoMk (imtrPreIsoOfIso P (pushId_uniq P F A hb Z).choose)
        (InducedWideCategory.Hom.ext (by
          -- ★`rw` は motive が壊れる(`pushHom` の**証明引数**が `Z.hom.hom` に依存する)。
          -- ★`congrArg` なら関数が明示なので依存が起きない。
          have h2 := congrArg (fun t => (pushId_uniq P F A hb Z).choose.hom ≫ t)
            (pushId_uniq P F A hb Z).choose_spec.1
          simp only [Iso.hom_inv_id_assoc] at h2
          -- ★`simpa` ではなく `exact`(既定透明度でないと 2 つの綴りが同一視されない、表 #1)
          exact h2)))
    (fun _ => imtrPreOver_hom_subsingleton P F _ _)

include F in
/-- ★★**`𝒪^×(A)^{imtr-pre}` は単位元を含む**(部分群性の 3 条件のうち 1 つ目)。 -/
theorem one_mem_otimesImtrPre (A : C) : (1 : End A) ∈ OTimesImtrPre P F A :=
  ⟨(OTimes P A).one_mem, isBaseIsomorphism_of_isIso P (𝟙 A),
   ⟨pushFunctorIdIso P F A (isBaseIsomorphism_of_isIso P (𝟙 A))⟩⟩

include F in
/-- ★★**合成の分解** —— `ε ≫ (φ ≫ ψ)` の**第2の分解**。

★左因子は co-angular base-isomorphism の**合成**(`F.coAngularComp` と
`isBaseIsomorphism_comp`)、右因子は isometric pre-step。 -/
theorem pushComp_uniq {A B E : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    (ψ : B ⟶ E) (hψ : IsBaseIsomorphism P ψ) (Z : Over (⟨A⟩ : ImtrPre P)) :
    ∃ δ : pushObj P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ)
            Z.hom.hom Z.hom.property.1 Z.hom.property.2
        ≅ pushObj P F ψ hψ
            (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
            (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
            (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2,
      pushHom P F ψ hψ
          (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
          (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
          (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
        = δ.inv ≫ pushHom P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ)
            Z.hom.hom Z.hom.property.1 Z.hom.property.2 ∧
      (pushMid P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2
        ≫ pushMid P F ψ hψ
            (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
            (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
            (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2)
        = pushMid P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ)
            Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ δ.hom := by
  set ε := Z.hom.hom with hε
  set hεi := Z.hom.property.1
  set hεs := Z.hom.property.2
  set h1 := (pushHom_spec P F φ hφ ε hεi hεs).1
  set h2 := (pushHom_spec P F φ hφ ε hεi hεs).2
  refine prop_1_9_i_uniq P F _ _ _ _
    (pushMid P F φ hφ ε hεi hεs ≫ pushMid P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2)
    (pushHom P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2) ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_
  · rw [← pushFac P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs,
      Category.assoc, ← pushFac P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2,
      ← Category.assoc, pushFac P F φ hφ ε hεi hεs, Category.assoc]
  · exact (pushMid_spec P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs).1
  · exact (pushMid_spec P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs).2
  · exact (pushHom_spec P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs).1
  · exact (pushHom_spec P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs).2
  · exact F.coAngularComp _ _ (pushMid_spec P F φ hφ ε hεi hεs).1
      (pushMid_spec P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2).1
  · exact isBaseIsomorphism_comp P (pushMid_spec P F φ hφ ε hεi hεs).2
      (pushMid_spec P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2).2
  · exact (pushHom_spec P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2).1
  · exact (pushHom_spec P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2).2

include F in
/-- ★★**`(φ ≫ ψ)_* ≅ φ_* ⋙ ψ_*`**。★自然性は hom の subsingleton 性から自動。 -/
noncomputable def pushFunctorCompIso {A B E : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    (ψ : B ⟶ E) (hψ : IsBaseIsomorphism P ψ) :
    pushFunctor P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ)
      ≅ pushFunctor P F φ hφ ⋙ pushFunctor P F ψ hψ :=
  NatIso.ofComponents
    (fun Z =>
      Over.isoMk (imtrPreIsoOfIso P (pushComp_uniq P F φ hφ ψ hψ Z).choose)
        (InducedWideCategory.Hom.ext (by
          have h2 := congrArg (fun t => (pushComp_uniq P F φ hφ ψ hψ Z).choose.hom ≫ t)
            (pushComp_uniq P F φ hφ ψ hψ Z).choose_spec.1
          simp only [Iso.hom_inv_id_assoc] at h2
          exact h2)))
    (fun _ => imtrPreOver_hom_subsingleton P F _ _)

include F in
/-- ★★**積で閉じる**(部分群性の 2 つ目)。

★**`End A` の積は `x * y = y ≫ x`** である(合成の向きに注意)。 -/
theorem mul_mem_otimesImtrPre (A : C) {u v : End A}
    (hu : u ∈ OTimesImtrPre P F A) (hv : v ∈ OTimesImtrPre P F A) :
    u * v ∈ OTimesImtrPre P F A := by
  obtain ⟨huo, hbu, ⟨eu⟩⟩ := hu
  obtain ⟨hvo, hbv, ⟨ev⟩⟩ := hv
  refine ⟨(OTimes P A).mul_mem huo hvo, isBaseIsomorphism_comp P hbv hbu, ⟨?_⟩⟩
  exact pushFunctorCompIso P F _ hbv _ hbu ≪≫
    Functor.isoWhiskerRight ev (pushFunctor P F (u : A ⟶ A) hbu) ≪≫
    Functor.leftUnitor _ ≪≫ eu

include F in
/-- ★★**逆元で閉じる**(部分群性の 3 つ目)。

★`u * w = 1`(すなわち `w ≫ u = 𝟙`)なる `w` が `𝒪^×(A)` にあれば、
`u ∈ ⟹ w ∈`。★**積と単位元の結果から出る。** -/
theorem inv_mem_otimesImtrPre (A : C) {u w : End A}
    (hu : u ∈ OTimesImtrPre P F A) (hwo : w ∈ OTimes P A) (h : u * w = 1) :
    w ∈ OTimesImtrPre P F A := by
  obtain ⟨-, hbu, ⟨eu⟩⟩ := hu
  -- `w` が base-isomorphism であること —— `𝒪^×` の元は同型
  haveI hwi : IsIso (w : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso (w : End A)).mp hwo.2
  have hbw : IsBaseIsomorphism P (w : A ⟶ A) := isBaseIsomorphism_of_isIso P _
  refine ⟨hwo, hbw, ⟨?_⟩⟩
  -- `w ≫ u = 𝟙` から `pushFunctor w ⋙ pushFunctor u ≅ 𝟭`
  -- ★`End A` の積は `x * y = y ≫ x` なので、`u * w = 1` は `w ≫ u = 𝟙` である
  have hcomp : ((w : A ⟶ A) ≫ (u : A ⟶ A)) = 𝟙 A := h
  have hEq : pushFunctor P F ((w : A ⟶ A) ≫ (u : A ⟶ A)) (isBaseIsomorphism_comp P hbw hbu)
      = pushFunctor P F (𝟙 A) (isBaseIsomorphism_of_isIso P (𝟙 A)) := by
    -- ★`rw` は motive が壊れる(証明引数が射に依存)。`congr` なら証明無関係が効く。
    congr 1
  have e1 : pushFunctor P F (w : A ⟶ A) hbw ⋙ pushFunctor P F (u : A ⟶ A) hbu
      ≅ 𝟭 (Over (⟨A⟩ : ImtrPre P)) :=
    (pushFunctorCompIso P F _ hbw _ hbu).symm ≪≫ eqToIso hEq ≪≫
      pushFunctorIdIso P F A (isBaseIsomorphism_of_isIso P (𝟙 A))
  exact (Functor.rightUnitor _).symm ≪≫ Functor.isoWhiskerLeft _ eu.symm ≪≫ e1

/-! ### ★★`u^{imtr-pre} = id` の**十分条件**

★**`u` が isometric pre-step と「同型を除いて可換」なら `u_* ≅ 𝟭`。**

★これは `Proposition 1.13, (iii)` で使う ——
`𝒞_A → 𝒞` の自己同型 `α` は、自然性から
`γ ≫ α_X = α_Y ≫ γ`(`γ` は isometric pre-step)を満たし、
`α_Y` は同型だから、まさにこの形になる。 -/

include F in
/-- ★`u` と `γ` が同型 `v` を挟んで可換なら、`γ ≫ u` の第2の分解が `(v, γ)` である。 -/
theorem pushComm_uniq {A : C} (u : A ⟶ A) (hb : IsBaseIsomorphism P u)
    (Z : Over (⟨A⟩ : ImtrPre P)) (v : Z.left.obj ⟶ Z.left.obj) [IsIso v]
    (hcomm : Z.hom.hom ≫ u = v ≫ Z.hom.hom) :
    ∃ δ : pushObj P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≅ Z.left.obj,
      Z.hom.hom = δ.inv ≫ pushHom P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ∧
      v = pushMid P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ δ.hom :=
  prop_1_9_i_uniq P F _ _ _ _ v Z.hom.hom
    (by rw [← pushFac P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2]; exact hcomm)
    (pushMid_spec P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushMid_spec P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushHom_spec P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushHom_spec P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (isCoAngular_of_isIso P v) (isBaseIsomorphism_of_isIso P v)
    Z.hom.property.1 Z.hom.property.2

include F in
/-- ★★**`u` が isometric pre-step と同型を除いて可換なら `u_* ≅ 𝟭`**。

★`pushFunctorIdIso` の一般化(`v = 𝟙` がその場合)。自然性は
`imtrPreOver_hom_subsingleton` から自動。 -/
noncomputable def pushFunctorIsoOfCommutes {A : C} (u : A ⟶ A) (hb : IsBaseIsomorphism P u)
    (comm : ∀ Z : Over (⟨A⟩ : ImtrPre P),
      Σ' v : Z.left.obj ⟶ Z.left.obj, PLift (IsIso v) ×' Z.hom.hom ≫ u = v ≫ Z.hom.hom) :
    pushFunctor P F u hb ≅ 𝟭 (Over (⟨A⟩ : ImtrPre P)) :=
  NatIso.ofComponents
    (fun Z =>
      haveI : IsIso (comm Z).1 := (comm Z).2.1.down
      Over.isoMk (imtrPreIsoOfIso P (pushComm_uniq P F u hb Z (comm Z).1 (comm Z).2.2).choose)
        (InducedWideCategory.Hom.ext (by
          have h2 := congrArg
            (fun t => (pushComm_uniq P F u hb Z (comm Z).1 (comm Z).2.2).choose.hom ≫ t)
            (pushComm_uniq P F u hb Z (comm Z).1 (comm Z).2.2).choose_spec.1
          simp only [Iso.hom_inv_id_assoc] at h2
          exact h2)))
    (fun _ => imtrPreOver_hom_subsingleton P F _ _)

include F in
/-- ★★**`𝒪^×(A)^{imtr-pre}` を部分モノイドとして包む**。

原文 (FrdI p.31):
> the subgroup of v ∈O×(A) for which vimtr-pre is the identity.

★**`𝒪^×(A)` 自体が `End A` の部分モノイドとして実装されている**ので、
ここも部分モノイドとして構成し、★**逆元の閉性は `inv_mem_otimesImtrPre` で別に述べる**。
`𝒪^×` の元はすべて可逆なので、この 3 条件が原文の「subgroup」の内容である。 -/
def OTimesImtrPreSubmonoid (A : C) : Submonoid (End A) where
  carrier := OTimesImtrPre P F A
  one_mem' := one_mem_otimesImtrPre P F A
  mul_mem' hx hy := mul_mem_otimesImtrPre P F A hx hy

include F in
/-- ★`𝒪^×(A)^{imtr-pre} ≤ 𝒪^×(A)` —— 原文の包含。 -/
theorem otimesImtrPreSubmonoid_le (A : C) :
    OTimesImtrPreSubmonoid P F A ≤ OTimes P A :=
  fun _ h => h.1

end ImtrPreFunctor

/-! ## ★第5段 —— (ii) の圏同値

原文 (FrdI p.31):
> C →A with φ : A →B. Moreover, if φ is a co-angular pre-step, then φ∗is an

★中身は原文 p.32–33 の 6 段である。まず `(Φ(v))⁻¹(Div v)` の計算規則を用意する。
-/

section PushEquiv

/-- ★**`(Φ(v))⁻¹(Div v) ∈ Φ(B)`** —— `Definition 1.3, (iii), (d)` のスライス側の不変量。

`coaPreOverFunctor` が対象に割り当てているものと**同じ式**である。 -/
noncomputable def divInv {X B : C} (v : X ⟶ B) (hv : IsBaseIsomorphism P v) :
    Φ.val (P.toElem.obj B).base :=
  haveI : IsIso (P.Base v) := hv
  Φ.map (inv (P.Base v)) (P.Div v)

/-- ★不変量は **isometric な射を左から合成しても変わらない**。

原文 (FrdI p.33) が「since `Div(β◦ψ) = Div(φ′◦α′)`, and `β`, `α′` are isometric」
と言って使うのがこれである。 -/
theorem divInv_comp_left {Y X B : C} (u : Y ⟶ X) (v : X ⟶ B)
    (hu : IsBaseIsomorphism P u) (hv : IsBaseIsomorphism P v) (hui : IsIsometric P u) :
    divInv P (u ≫ v) (isBaseIsomorphism_comp P hu hv) = divInv P v hv := by
  haveI hiu : IsIso (P.Base u) := hu
  haveI hiv : IsIso (P.Base v) := hv
  haveI : IsIso (P.Base (u ≫ v)) := isBaseIsomorphism_comp P hu hv
  show Φ.map (inv (P.Base (u ≫ v))) (P.Div (u ≫ v)) = Φ.map (inv (P.Base v)) (P.Div v)
  have hd : P.Div (u ≫ v) = Φ.map (P.Base u) (P.Div v) := by
    rw [P.Div_comp, show P.Div u = 0 from hui, smul_zero, add_zero]
  have hinv : inv (P.Base (u ≫ v)) = inv (P.Base v) ≫ inv (P.Base u) := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [P.Base_comp]
    simp
  rw [hd, hinv, Φ.map_comp (inv (P.Base u)) (inv (P.Base v)),
    ← Φ.map_comp (P.Base u) (inv (P.Base u)), IsIso.inv_hom_id, Φ.map_id]

/-- ★不変量に **isometric な pre-step を右から合成する**と、`Φ` で引き戻される。 -/
theorem divInv_comp_right {Y X B : C} (u : Y ⟶ X) (v : X ⟶ B)
    (hu : IsBaseIsomorphism P u) (hv : IsBaseIsomorphism P v) (hvi : IsIsometric P v)
    (hvl : IsLinear P v) :
    divInv P (u ≫ v) (isBaseIsomorphism_comp P hu hv) =
      haveI : IsIso (P.Base v) := hv
      Φ.map (inv (P.Base v)) (divInv P u hu) := by
  haveI hiu : IsIso (P.Base u) := hu
  haveI hiv : IsIso (P.Base v) := hv
  haveI : IsIso (P.Base (u ≫ v)) := isBaseIsomorphism_comp P hu hv
  show Φ.map (inv (P.Base (u ≫ v))) (P.Div (u ≫ v))
    = Φ.map (inv (P.Base v)) (Φ.map (inv (P.Base u)) (P.Div u))
  have hd : P.Div (u ≫ v) = P.Div u := by
    rw [P.Div_comp, show P.Div v = 0 from hvi, map_zero, zero_add,
      show P.degFr v = 1 from hvl]
    simp
  have hinv : inv (P.Base (u ≫ v)) = inv (P.Base v) ≫ inv (P.Base u) := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [P.Base_comp]
    simp
  rw [hd, hinv, Φ.map_comp (inv (P.Base u)) (inv (P.Base v))]

/-- ★`Φ` が divisorial なので、`Order(Φ(B))` の同型は**等式**に落ちる。

★`mle_antisymm` の 2 仮定(integral / sharp)は**どちらも `divisorial` に含まれる**。 -/
theorem eq_of_orderCat_iso {B : C} {a b : Φ.val (P.toElem.obj B).base}
    (e : toOrderCat a ≅ toOrderCat b) : a = b :=
  mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2 (leOfHom e.hom) (leOfHom e.inv)

/-- ★射が等しければ不変量も等しい(証明の取り方には依らない)。 -/
theorem divInv_congr {X B : C} {v v' : X ⟶ B} (h : v = v')
    (hv : IsBaseIsomorphism P v) (hv' : IsBaseIsomorphism P v') :
    divInv P v hv = divInv P v' hv' := by
  subst h; rfl

end PushEquiv

section PushEquiv2

variable [(coaPreProp P).IsMultiplicative]

/-- ★★**同じ不変量を持つ2つの co-angular pre-step は同型で移り合う**(スライス側)。

★`coaPre_iso_of_div_eq`(コスライス側、`Proposition 1.8` で作った)の**スライス版**。
`Definition 1.3, (iii), (d)` の**第2**の圏同値の充満性と忠実性を使う。 -/
theorem coaPreOver_iso_of_divInv_eq
    (hover : ∀ X : C, (coaPreOverFunctor P X).IsEquivalence)
    {B : C} (Z W : Over (⟨B⟩ : WideSubcategory (coaPreProp P)))
    (h : divInv P Z.hom.hom Z.hom.property.2.2 = divInv P W.hom.hom W.hom.property.2.2) :
    Nonempty (Z ≅ W) := by
  haveI := (hover B).full
  haveI := (hover B).faithful
  have hobj : (coaPreOverFunctor P B).obj Z = (coaPreOverFunctor P B).obj W :=
    congrArg (fun x => op (toOrderCat x)) h
  exact ⟨(coaPreOverFunctor P B).preimageIso (eqToIso hobj)⟩

/-- ★★**原文 p.32–33 の構成**。

`φ : A ⟶ B` が co-angular pre-step、`β : D ⟶ B` が isometric pre-step のとき、
co-angular pre-step `ψ : Cc ⟶ D` と isometric pre-step `α : Cc ⟶ A` で
`ψ ≫ β = α ≫ φ`(原文の `β ◦ ψ = φ ◦ α`)となるものが取れる。

★★これが (ii) の**本質的全射性**であり、原文が言うとおり
**`φ` を `ψ` に取り替えれば充満性**にもなる。

★使う入力: (iii)(d) の**第2**の圏同値(本質的全射性・充満性・忠実性)、
`Definition 1.3, (v), (c)`(`preStepFactor'`)、そして `Φ` の **divisorial 性**
(`Order(Φ(B))` の反対称性)。 -/
theorem push_square (F : FrobenioidCore P)
    (hover : ∀ X : C, (coaPreOverFunctor P X).IsEquivalence)
    {A B Dd : C} (φ : A ⟶ B) (hφc : IsCoAngular P φ) (hφs : IsPreStep P φ)
    (β : Dd ⟶ B) (hβi : IsIsometric P β) (hβs : IsPreStep P β) :
    ∃ (Cc : C) (ψ : Cc ⟶ Dd) (α : Cc ⟶ A),
      ψ ≫ β = α ≫ φ ∧ (IsCoAngular P ψ ∧ IsPreStep P ψ) ∧
        (IsIsometric P α ∧ IsPreStep P α) := by
  haveI := (hover Dd).essSurj
  haveI hbβ : IsIso (P.Base β) := hβs.2
  -- 1. 目標の不変量 `e ∈ Φ(D)`
  obtain ⟨e, he⟩ : ∃ e : Φ.val (P.toElem.obj Dd).base,
      e = Φ.map (P.Base β) (divInv P φ hφs.2) := ⟨_, rfl⟩
  -- 2. (iii)(d) の**本質的全射性**で co-angular pre-step `ψ` を取る
  obtain ⟨Z, hZ⟩ : ∃ Z : Over (⟨Dd⟩ : WideSubcategory (coaPreProp P)),
      Z = (coaPreOverFunctor P Dd).objPreimage (op (toOrderCat e)) := ⟨_, rfl⟩
  have hiZ : (coaPreOverFunctor P Dd).obj Z ≅ op (toOrderCat e) := by
    rw [hZ]; exact (coaPreOverFunctor P Dd).objObjPreimageIso _
  have hobj : (coaPreOverFunctor P Dd).obj Z
      = op (toOrderCat (divInv P Z.hom.hom Z.hom.property.2.2)) := rfl
  rw [hobj] at hiZ
  have hdivψ : divInv P Z.hom.hom Z.hom.property.2.2 = e :=
    (eq_of_orderCat_iso P hiZ.unop).symm
  -- 3. `ψ ≫ β` の不変量は `φ` のそれに等しい
  have hψβs : IsPreStep P (Z.hom.hom ≫ β) := IsPreStep.comp P Z.hom.property.2 hβs
  have hkey : divInv P (Z.hom.hom ≫ β)
      (isBaseIsomorphism_comp P Z.hom.property.2.2 hβs.2) = divInv P φ hφs.2 := by
    rw [divInv_comp_right P Z.hom.hom β Z.hom.property.2.2 hβs.2 hβi hβs.1, hdivψ, he,
      ← Φ.map_comp (P.Base β) (inv (P.Base β)), IsIso.inv_hom_id, Φ.map_id]
  -- 4. `Definition 1.3, (v), (c)` で `ψ ≫ β = α' ≫ φ'`
  obtain ⟨Xx, α', φ', hfac, hα'i, hα's, hφ'c, hφ's⟩ := F.preStepFactor' (Z.hom.hom ≫ β) hψβs
  -- 5. `φ'` の不変量も `φ` のそれに等しい
  have hdiv' : divInv P φ' hφ's.2 = divInv P φ hφs.2 := by
    rw [← divInv_comp_left P α' φ' hα's.2 hφ's.2 hα'i,
      divInv_congr P hfac.symm (isBaseIsomorphism_comp P hα's.2 hφ's.2)
        (isBaseIsomorphism_comp P Z.hom.property.2.2 hβs.2)]
    exact hkey
  -- 6. (iii)(d) の**充満性・忠実性**で同型 `γ : Xx ≅ A`
  obtain ⟨iso⟩ := coaPreOver_iso_of_divInv_eq P hover
    (Over.mk (⟨φ', hφ'c, hφ's⟩ :
      (⟨Xx⟩ : WideSubcategory (coaPreProp P)) ⟶ (⟨B⟩ : WideSubcategory (coaPreProp P))))
    (Over.mk (⟨φ, hφc, hφs⟩ :
      (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ (⟨B⟩ : WideSubcategory (coaPreProp P))))
    hdiv'
  obtain ⟨γ, hγ⟩ : ∃ γ : Xx ⟶ A, γ = (Over.Hom.left iso.hom).hom := ⟨_, rfl⟩
  obtain ⟨γ', hγ'⟩ : ∃ γ' : A ⟶ Xx, γ' = (Over.Hom.left iso.inv).hom := ⟨_, rfl⟩
  haveI hγiso : IsIso γ := by
    refine ⟨γ', ?_, ?_⟩
    · rw [hγ, hγ']
      exact congrArg (fun t => InducedWideCategory.Hom.hom (Over.Hom.left t)) iso.hom_inv_id
    · rw [hγ, hγ']
      exact congrArg (fun t => InducedWideCategory.Hom.hom (Over.Hom.left t)) iso.inv_hom_id
  have hγφ : γ ≫ φ = φ' := by
    rw [hγ]
    exact congrArg InducedWideCategory.Hom.hom (Over.w iso.hom)
  -- 7. `α := α' ≫ γ`
  refine ⟨Z.left.obj, Z.hom.hom, α' ≫ γ, ?_, Z.hom.property, ?_, ?_⟩
  · rw [Category.assoc, hγφ, hfac]
  · exact IsIsometric.comp P hα'i (isIsometric_of_isIso P γ)
  · exact IsPreStep.comp P hα's (isPreStep_of_isIso P γ)

/-- 補助: スライスの同型から `𝒞` の同型を取り出す。 -/
theorem coaPreOver_hom_of_iso {B : C} {Z W : Over (⟨B⟩ : WideSubcategory (coaPreProp P))}
    (iso : Z ≅ W) : ∃ γ : Z.left.obj ⟶ W.left.obj, IsIso γ ∧ γ ≫ W.hom.hom = Z.hom.hom := by
  refine ⟨(Over.Hom.left iso.hom).hom, ⟨(Over.Hom.left iso.inv).hom, ?_, ?_⟩, ?_⟩
  · exact congrArg (fun t => InducedWideCategory.Hom.hom (Over.Hom.left t)) iso.hom_inv_id
  · exact congrArg (fun t => InducedWideCategory.Hom.hom (Over.Hom.left t)) iso.inv_hom_id
  · exact congrArg InducedWideCategory.Hom.hom (Over.w iso.hom)

/-- ★(i) の分解の左因子は **pre-step でもある**(`Proposition 1.7, (v)`)。 -/
theorem pushMid_isPreStep (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) (hφl : IsLinear P φ) {Cc : C} (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    IsPreStep P (pushMid P F φ hφ ε hεi hεs) := by
  refine ⟨?_, (pushMid_spec P F φ hφ ε hεi hεs).2⟩
  have hlin : IsLinear P (ε ≫ φ) := IsLinear.comp P hεs.1 hφl
  rw [pushFac P F φ hφ ε hεi hεs] at hlin
  exact (prop_1_7_v_linear P _ _ hlin).1

/-- ★★★**`Proposition 1.9, (ii)` の圏同値** —— `φ` が co-angular pre-step なら
`φ_* : 𝒞^imtr-pre_A ⥤ 𝒞^imtr-pre_B` は圏同値。

★**忠実性は無料**(`𝒞^imtr-pre_A` が thin)。
★**本質的全射性**は `push_square` そのもの。
★**充満性**は `push_square` を **`φ` の代わりに `pushMid W`、`β` の代わりに `h`** に当てる ——
原文の「by possibly replacing φ by ψ」がこれである。 -/
theorem pushFunctor_isEquivalence (F : FrobenioidCore P)
    (hover : ∀ X : C, (coaPreOverFunctor P X).IsEquivalence)
    {A B : C} (φ : A ⟶ B) (hφc : IsCoAngular P φ) (hφs : IsPreStep P φ) :
    (pushFunctor P F φ hφs.2).IsEquivalence := by
  haveI hmono : Mono φ := F.preStepMono φ hφs
  haveI hfaith : (pushFunctor P F φ hφs.2).Faithful :=
    ⟨fun _ => imtrPreOver_hom_subsingleton P F _ _⟩
  -- ★★充満性
  haveI hfull : (pushFunctor P F φ hφs.2).Full := by
    constructor
    intro Z W h
    -- 記号
    obtain ⟨hh, hhe⟩ : ∃ hh : pushObj P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ⟶
        pushObj P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2,
        hh = (Over.Hom.left h).hom := ⟨_, rfl⟩
    have hhp : IsIsometric P hh ∧ IsPreStep P hh := by
      rw [hhe]; exact (Over.Hom.left h).property
    have hhw : hh ≫ pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 =
        pushHom P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 := by
      rw [hhe]; exact congrArg InducedWideCategory.Hom.hom (Over.w h)
    have hmW := pushMid_spec P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2
    have hmZ := pushMid_spec P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2
    have hmWs := pushMid_isPreStep P F φ hφs.2 hφs.1 W.hom.hom W.hom.property.1 W.hom.property.2
    have hmZs := pushMid_isPreStep P F φ hφs.2 hφs.1 Z.hom.hom Z.hom.property.1 Z.hom.property.2
    have hHW := pushHom_spec P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2
    -- ★`push_square` を `φ := pushMid W`, `β := hh` に当てる
    obtain ⟨Cc', ψ', a, heq', hψ', ha⟩ := push_square P F hover
      (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) hmW.1 hmWs
      hh hhp.1 hhp.2
    -- ★不変量の一致
    have d1 : divInv P (ψ' ≫ hh) (isBaseIsomorphism_comp P hψ'.2.2 hhp.2.2)
        = divInv P (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) hmW.2 := by
      rw [divInv_congr P heq' (isBaseIsomorphism_comp P hψ'.2.2 hhp.2.2)
        (isBaseIsomorphism_comp P ha.2.2 hmW.2)]
      exact divInv_comp_left P a _ ha.2.2 hmW.2 ha.1
    have d4 : divInv P (pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ hh)
        (isBaseIsomorphism_comp P hmZ.2 hhp.2.2)
        = divInv P (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) hmW.2 := by
      haveI : IsIso (P.Base (pushHom P F φ hφs.2 W.hom.hom
        W.hom.property.1 W.hom.property.2)) := hHW.2.2
      refine Φ.map_injective (inv (P.Base (pushHom P F φ hφs.2 W.hom.hom
        W.hom.property.1 W.hom.property.2))) ?_
      have e1 := divInv_comp_right P
        (pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ hh)
        (pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2)
        (isBaseIsomorphism_comp P hmZ.2 hhp.2.2) hHW.2.2 hHW.1 hHW.2.1
      have e2 := divInv_comp_right P
        (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2)
        (pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2)
        hmW.2 hHW.2.2 hHW.1 hHW.2.1
      rw [← e1, ← e2]
      have f1 : (pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ hh) ≫
          pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 = Z.hom.hom ≫ φ := by
        rw [Category.assoc, hhw, ← pushFac P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2]
      have f2 : pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 ≫
          pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 = W.hom.hom ≫ φ :=
        (pushFac P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2).symm
      rw [divInv_congr P f1 _ (isBaseIsomorphism_comp P Z.hom.property.2.2 hφs.2),
        divInv_congr P f2 _ (isBaseIsomorphism_comp P W.hom.property.2.2 hφs.2),
        divInv_comp_left P Z.hom.hom φ Z.hom.property.2.2 hφs.2 Z.hom.property.1,
        divInv_comp_left P W.hom.hom φ W.hom.property.2.2 hφs.2 W.hom.property.1]
    have dkey : divInv P (pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2) hmZ.2
        = divInv P ψ' hψ'.2.2 := by
      haveI : IsIso (P.Base hh) := hhp.2.2
      refine Φ.map_injective (inv (P.Base hh)) ?_
      rw [← divInv_comp_right P _ hh hmZ.2 hhp.2.2 hhp.1 hhp.2.1,
        ← divInv_comp_right P _ hh hψ'.2.2 hhp.2.2 hhp.1 hhp.2.1, d4, ← d1]
    -- ★同型 `ee : Z.left.obj ≅ Cc'`
    obtain ⟨iso⟩ := coaPreOver_iso_of_divInv_eq P hover
      (Over.mk (⟨pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2, hmZ.1, hmZs⟩ :
        (⟨Z.left.obj⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨_⟩))
      (Over.mk (⟨ψ', hψ'.1, hψ'.2⟩ : (⟨Cc'⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨_⟩))
      dkey
    obtain ⟨ee, hee, heew0⟩ := coaPreOver_hom_of_iso P iso
    have heew : ee ≫ ψ' = pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 := heew0
    -- ★求める射
    refine ⟨Over.homMk (⟨ee ≫ a, IsIsometric.comp P (isIsometric_of_isIso P ee) ha.1,
        IsPreStep.comp P (isPreStep_of_isIso P ee) ha.2⟩ : Z.left ⟶ W.left) ?_,
      imtrPreOver_hom_subsingleton P F _ _⟩
    refine InducedWideCategory.Hom.ext ?_
    show (ee ≫ a) ≫ W.hom.hom = Z.hom.hom
    refine (cancel_mono φ).mp ?_
    calc ((ee ≫ a) ≫ W.hom.hom) ≫ φ = (ee ≫ a) ≫ (W.hom.hom ≫ φ) := Category.assoc _ _ _
      _ = (ee ≫ a) ≫ (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) := by
          rw [← pushFac P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2]
      _ = ee ≫ ((a ≫ pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) := by
          simp only [Category.assoc]
      _ = ee ≫ ((ψ' ≫ hh) ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) :=
          congrArg (fun t => ee ≫ (t ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2)) heq'.symm
      _ = (ee ≫ ψ') ≫ (hh ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) := by
          simp only [Category.assoc]
      _ = pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫
            pushHom P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 := by
          rw [heew, hhw]
          rfl
      _ = Z.hom.hom ≫ φ :=
          (pushFac P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2).symm
  -- ★★本質的全射性
  haveI hess : (pushFunctor P F φ hφs.2).EssSurj := by
    constructor
    intro Y
    obtain ⟨Cc, ψ, α, heq, hψ, hα⟩ := push_square P F hover φ hφc hφs
      Y.hom.hom Y.hom.property.1 Y.hom.property.2
    obtain ⟨δ, hδ1, -⟩ := prop_1_9_i_uniq P F
      (pushObj P F φ hφs.2 α hα.1 hα.2) Y.left.obj
      (pushMid P F φ hφs.2 α hα.1 hα.2) (pushHom P F φ hφs.2 α hα.1 hα.2)
      ψ Y.hom.hom
      ((pushFac P F φ hφs.2 α hα.1 hα.2).symm.trans heq.symm)
      (pushMid_spec P F φ hφs.2 α hα.1 hα.2).1 (pushMid_spec P F φ hφs.2 α hα.1 hα.2).2
      (pushHom_spec P F φ hφs.2 α hα.1 hα.2).1 (pushHom_spec P F φ hφs.2 α hα.1 hα.2).2
      hψ.1 hψ.2.2 Y.hom.property.1 Y.hom.property.2
    refine ⟨Over.mk (⟨α, hα⟩ : (⟨Cc⟩ : ImtrPre P) ⟶ (⟨A⟩ : ImtrPre P)), ⟨?_⟩⟩
    refine Over.isoMk ⟨⟨δ.hom, isIsometric_of_isIso P δ.hom, isPreStep_of_isIso P δ.hom⟩,
      ⟨δ.inv, isIsometric_of_isIso P δ.inv, isPreStep_of_isIso P δ.inv⟩,
      InducedWideCategory.Hom.ext δ.hom_inv_id,
      InducedWideCategory.Hom.ext δ.inv_hom_id⟩ ?_
    refine InducedWideCategory.Hom.ext ?_
    show δ.hom ≫ Y.hom.hom = pushHom P F φ hφs.2 α hα.1 hα.2
    rw [hδ1, ← Category.assoc, δ.hom_inv_id, Category.id_comp]
  exact ⟨hfaith, hfull, hess⟩

end PushEquiv2

/-! ## ★負の対照 —— (iv) の `co-angular` は落とせない

★(iv) は「co-angular **かつ** linear なら isotropic 性が両向きに移る」と言う。
`linear` だけでは `⇐` が壊れることを、`Proposition 1.4` で作った模型 `cP` で示す。

★`cP` は `Vee` 上の定数関手による pre-Frobenioid で、
**`vLeftTop : left ⟶ top` は linear だが co-angular でなく**、
**`top` は isotropic だが `left` は isotropic でない**。 -/

/-- `cP` の `vLeftTop` は **linear**。 -/
theorem cP_linear_vLeftTop : IsLinear cP vLeftTop := rfl

/-- `cP` の `top` は **isotropic** —— `top` から出る射は `𝟙` しかない。 -/
theorem cP_isotropic_top : IsIsotropic cP Vee.top := by
  intro Dd f _ _
  have hd : Dd = Vee.top := by
    rcases leOfHom f with h | h
    · exact h.symm
    · exact h
  subst hd
  haveI := Preorder.subsingleton_hom Vee.top Vee.top
  exact ⟨𝟙 _, Subsingleton.elim _ _, Subsingleton.elim _ _⟩

/-- ★★**(iv) の `co-angular` は落とせない** —— `vLeftTop` は linear で
終域 `top` は isotropic だが、始域 `left` は **isotropic でない**。

★`Proposition 1.4, (i)` の負の対照(`cP_not_coAngular`)がそのまま
`Proposition 1.9, (iv)` の負の対照になる。**同じ模型が2つの命題で効いた。** -/
theorem cP_prop_1_9_iv_coAngular_is_load_bearing :
    IsLinear cP vLeftTop ∧ IsIsotropic cP Vee.top ∧
      ¬ IsCoAngular cP vLeftTop ∧ ¬ IsIsotropic cP Vee.left :=
  ⟨cP_linear_vLeftTop, cP_isotropic_top, cP_not_coAngular, cP_not_isotropic_left⟩

/-! ### ★出典の紐付け(`.src`)

★`.src` は「その原典項目を**完全に**実装した」という主張である。
`Proposition 1.9` は (i)–(vii) がすべて揃ったので、ここで初めて付ける。
各条の内訳:

* **(i)**: `prop_1_9_i_factor`(分解)＋ `prop_1_9_i_uniq`(一意性)＋
  `prop_1_9_i_isometric_iff` / `prop_1_9_i_coAngular_iff` / `prop_1_9_i_pullBack_iff`
* **(ii)**: `pushFunctor`(関手)＋ `pushFunctor_isEquivalence`(圏同値)＋
  `OTimesImtrPre`(`𝒪^×(A)^imtr-pre`)
* **(iii)**: `pullFunctor`(関手)＋ `prop_1_9_iii_uniq`(同型を除く一意性)
* **(iv)**: `prop_1_9_iv`
* **(v)**: `istr_frobenioid`(`𝒞^istr` が Frobenioid)＋ `isotropificationAdj`(左随伴)＋
  `istr_compat_degFr` / `_Base` / `_Div`(`𝒞 → 𝔽_Φ` が経由すること)＋
  `isotropificationRestrictIso`(制限が恒等と同型)＋ 保存 11 種
  (`isotropification_frobType` / `_degFr` / `_preStep_iff` / `_pullBack` / `_baseIso_iff` /
  `_baseFSM` / `_baseIdentity` / `_divIdentity` / `_isometric_iff` / `istr_coAngular` /
  `_lbInvertible`)
* **(vi)**: `prop_1_9_vi`
* **(vii)**: `prop_1_9_vii`
-/

def prop_1_9_i_factor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 31, item := "Proposition 1.9, (i)",
    sectionId := "frdi-prop-1-9-i" }

def pushFunctor_isEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 31, item := "Proposition 1.9, (ii)",
    sectionId := "frdi-prop-1-9-ii" }

def pullFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 31, item := "Proposition 1.9, (iii)",
    sectionId := "frdi-prop-1-9-iii" }

def prop_1_9_iv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 31, item := "Proposition 1.9, (iv)",
    sectionId := "frdi-prop-1-9-iv" }

def istr_frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 32, item := "Proposition 1.9, (v)",
    sectionId := "frdi-prop-1-9-v" }

def prop_1_9_vi.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 32, item := "Proposition 1.9, (vi)",
    sectionId := "frdi-prop-1-9-vi" }

def prop_1_9_vii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 32, item := "Proposition 1.9, (vii)",
    sectionId := "frdi-prop-1-9-vii" }

/-! ## ★★命題全体の `.src`(2026-08-16)

★上の `.src` はすべて**条つき**で、locator の記録である(器具の数には入らない)。
★下の `prop_1_9.src` が**命題全体を完全に実装したという主張**である。

★★**文脈を持たない検証役の監査(2026-08-16)で欠落が 2 件見つかり、どちらも同日に埋めた。**

原文 (FrdI p.32):
> forms a left adjoint to the inclusion functor Cistr →C, through which the functor

★★**(v) の 2 件**:
- 「through which the functor `𝒞 → 𝔽_Φ` factors」→ `isotropificationFactorIso`
  (`P.toElem ≅ isotropification P F ⋙ (istrPre P F).toElem`)。
  ★中身は「**isotropic hull は `𝔽_Φ` では同型になる**」(`toElem_map_hullMap_isIso`) ——
  isometric ＋ pre-step から `div = 0`・`deg = 1`・base 同型が揃うため。
- 保存 11 クラスのうち co-angular → `isotropification_coAngular`(仮定不要。原文より強い)。

★★**(ii) の 1 件**:

原文 (FrdI p.31):
> the subgroup of v ∈O×(A) for which vimtr-pre is the identity.

★**原文は `𝒪^×(A)^{imtr-pre}` が `𝒪^×(A)` の部分群であると述べている**が、
監査以前は `Set (End A)` と包含補題しか無かった。3 条件を埋めた:
- 単位元: `one_mem_otimesImtrPre`(`ε` 自身が isometric pre-step なので `ε = 𝟙 ≫ ε` が第2の分解)
- 積: `mul_mem_otimesImtrPre`(`ε ≫ (φ ≫ ψ)` の第2の分解 —— 左因子は
  co-angular base-isomorphism の**合成**、右因子は isometric pre-step)
- 逆元: `inv_mem_otimesImtrPre`(積と単位元から)
- 包装: `OTimesImtrPreSubmonoid`, `otimesImtrPreSubmonoid_le`

★★**共通の鍵**: `Over (⟨A⟩ : ImtrPre P)` の hom は subsingleton
(`imtrPreOver_hom_subsingleton`)。★**したがって自然性は自動**で、
対象ごとの同型さえ作れば関手の同型が得られる。

## ★条ごとの照合表

| 条 | 主張 | 宣言 |
|---|---|---|
| (i) | 分解の存在 | `prop_1_9_i_factor` |
| (i) | 一意性(`(α◦γ, γ⁻¹◦β)` を除く) | `prop_1_9_i_uniq` |
| (i) | φ isometric ⟺ β が Frobenius 型 | `prop_1_9_i_isometric_iff` |
| (i) | φ co-angular ⟺ α が同型 | `prop_1_9_i_coAngular_iff` |
| (i) | φ pull-back ⟺ φ が同型 | `prop_1_9_i_pullBack_iff` |
| (ii) | 関手 `φ_*` | `pushFunctor` |
| (ii) | **Moreover**: co-angular pre-step なら圏同値 | `pushFunctor_isEquivalence` |
| (ii) | `𝒪^×(A)^{imtr-pre}` と包含 | `OTimesImtrPre`, `otimesImtrPre_subset` |
| (ii) | ★**部分群であること**(監査で発見) | `OTimesImtrPreSubmonoid` ほか 3 本 |
| (iii) | 関手 `φ^*` | `pullFunctor` |
| (iii) | 四角形と、同型を除く一意性 | `pullFac`, `prop_1_9_iii_uniq`, `prop_1_9_iii_lift` |
| (iv) | iff | `prop_1_9_iv` |
| (v) | `𝒞^istr` が Frobenioid | `istr_frobenioid` |
| (v) | isotropification が左随伴 | `isotropificationAdj` |
| (v) | ★**`𝒞 → 𝔽_Φ` が経由する**(監査で発見) | `isotropificationFactorIso` |
| (v) | 制限が恒等と同型 | `isotropificationRestrictIso` |
| (v) | 保存 11 クラス(★`isotropification_coAngular` は監査で追加) | `isotropification_*` 11 本 |
| (v) | moreover: 包含関手との互換性 | `istr_compat_*`, `istr_frobType_iff`, `istr_isPullBack_*` |
| (vi) | iff ＋ moreover | `prop_1_9_vi` |
| (vii) | iff | `prop_1_9_vii` |

★**★印の 3 行が、監査以前には存在しなかった主張である。**
-/

def prop_1_9.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 30, item := "Proposition 1.9",
    sectionId := "frdi-prop-1-9-i" }

end ABC3.Found.FrdI
