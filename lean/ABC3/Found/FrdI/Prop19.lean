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

/-- ★★**`𝒞^istr` の pre-Frobenioid 構造** —— 原文の
「equipped with the **restriction** to `𝒞` of the given functor `𝒞 → 𝔽_Φ`」。

★★**4フィールドのうち3つが `P` のものそのまま**である。
`totEpiC` だけが1行の議論を要した。**これが「移送」の意味である。** -/
def istrPre : PreFrobenioid (Istr P) Φ where
  toElem := (isotropicProp P).ι ⋙ P.toElem
  divisorial := P.divisorial
  totEpiC := istr_totEpi P
  totEpiD := P.totEpiD

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
    (istrPre P).degFr g = P.degFr g.hom := rfl

theorem istr_compat_Base {X Y : Istr P} (g : X ⟶ Y) :
    (istrPre P).Base g = P.Base g.hom := rfl

theorem istr_compat_Div {X Y : Istr P} (g : X ⟶ Y) :
    (istrPre P).Div g = P.Div g.hom := rfl

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

theorem istr_isotropic (X : Istr P) : IsIsotropic (istrPre P) X := by
  intro Dd φ hi hs
  haveI : IsIso φ.hom := X.property Dd.obj φ.hom hi hs
  exact ⟨InducedCategory.homMk (inv φ.hom), InducedCategory.hom_ext (by simp),
    InducedCategory.hom_ext (by simp)⟩

/-- ★★**`𝒞^istr` のすべての射は co-angular**(`Proposition 1.4, (i)`)。

★原文が保存リストで「co-angular morphisms [cf. Proposition 1.4, (i)]」と
括弧書きするのはこれ —— **`𝒞^istr` では co-angular 性が自明になる**ので、
「保存する」は言うまでもない。 -/
theorem istr_coAngular {X Y : Istr P} (g : X ⟶ Y) : IsCoAngular (istrPre P) g :=
  prop_1_4_i (istrPre P) g (fun Z _ => istr_isotropic P Z)

/-! ### ★★21 条の「移送」—— まず辞書と、移送しやすい条から

★`Definition 1.3` の各条を `istrPre P` について示す。原文は
「[from the fact that `𝒞` is a Frobenioid!]」の一言だが、
**Lean では `Istr P` と `C` の間の辞書を1本ずつ引く**必要がある。
どこまでが本当に「移送」かを条ごとに測る。 -/

include F in
/-- ★辞書: `Istr P` の Frobenius 型は `C` のそれ。

★co-angular は**両側で自動**なので、実質は isometric と base-isomorphism の移送。 -/
theorem istr_frobType_iff {X Y : Istr P} (g : X ⟶ Y) :
    IsFrobeniusType (istrPre P) g ↔ IsFrobeniusType P g.hom :=
  ⟨fun h => ⟨⟨isCoAngular_of_isotropic_dom P F X.property g.hom, h.1.2⟩, h.2⟩,
   fun h => ⟨⟨istr_coAngular P g, h.1.2⟩, h.2⟩⟩

include F in
/-- **(iii)(a)** の移送 —— `𝒞^istr` では co-angular が自動なので**自明**。 -/
theorem istr_coAngularComp {X Y Z : Istr P} (ψ : X ⟶ Y) (φ : Y ⟶ Z) :
    IsCoAngular (istrPre P) ψ → IsCoAngular (istrPre P) φ →
      IsCoAngular (istrPre P) (ψ ≫ φ) :=
  fun _ _ => istr_coAngular P _

include F in
/-- **(iii)(b)** の移送 —— 同上。 -/
theorem istr_coAngularOfPreStep {X Y : Istr P} (α : X ⟶ Y) :
    IsCoAngular (istrPre P) α → IsPreStep (istrPre P) α →
      ∀ φ : X ⟶ Y, IsCoAngular (istrPre P) φ :=
  fun _ _ φ => istr_coAngular P φ

include F in
/-- **(vii)(b)** の移送 —— `𝒞^istr` の対象はすべて isotropic なので**自明**。 -/
theorem istr_isotropicClosed {X Y : Istr P} (_φ : X ⟶ Y) :
    IsIsotropic (istrPre P) X → IsIsotropic (istrPre P) Y :=
  fun _ => istr_isotropic P Y

include F in
/-- **(vii)(a)** の移送 —— `X` 自身が isotropic なので `𝟙_X` が isotropic hull。 -/
theorem istr_isotropicHullExists (X : Istr P) :
    ∃ (Y : Istr P) (φ : X ⟶ Y), IsIsotropicHull (istrPre P) φ :=
  ⟨X, 𝟙 X, (istrPre P).Div_id X, isPreStep_id _ X, istr_isotropic P X,
    fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm, fun β hβ => by
      have hg : γ = β := by simpa using hβ
      exact hg.symm⟩⟩

include F in
/-- **(v)(a)** の移送 —— `C` の mono 性が充満部分圏へそのまま降りる。 -/
theorem istr_preStepMono {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P) φ) :
    Mono φ := by
  haveI : Mono φ.hom := F.preStepMono φ.hom hφ
  refine ⟨fun {Z} g h hgh => ?_⟩
  refine InducedCategory.hom_ext ?_
  exact (cancel_mono φ.hom).mp (congrArg InducedCategory.Hom.hom hgh)

include F in
/-- **(ii)** の本質的一意性の移送 —— `C` で得た同型を充満部分圏へ持ち上げるだけ。 -/
theorem istr_frobDegUniq (X Y Z : Istr P) (φ : X ⟶ Y) (ψ : X ⟶ Z)
    (hφ : IsFrobeniusType (istrPre P) φ) (hψ : IsFrobeniusType (istrPre P) ψ)
    (hd : (istrPre P).degFr φ = (istrPre P).degFr ψ) :
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
    IsPullBack (istrPre P) g := by
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
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (istrPre P) γ ∧
        IsPreStep (istrPre P) β ∧ IsPullBack (istrPre P) α := by
  obtain ⟨Z₀, W₀, γ₀, β₀, α₀, heq, hγ, hβ, hα⟩ := F.arbFactor φ.hom
  have hZ : IsIsotropic P Z₀ := F.isotropicClosed γ₀ X.property
  have hW : IsIsotropic P W₀ := F.isotropicClosed β₀ hZ
  refine ⟨⟨Z₀, hZ⟩, ⟨W₀, hW⟩, InducedCategory.homMk γ₀, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, ?_, hβ,
    istr_isPullBack_of P _ hα⟩
  exact (istr_frobType_iff P F (X := X) (Y := ⟨Z₀, hZ⟩)
    (InducedCategory.homMk γ₀)).mpr hγ

include F in
/-- **(ii)** の移送 —— 各次数の Frobenius 型射。中間対象は自動で isotropic。 -/
theorem istr_frobDegSurj (X : Istr P) (n : ℕ+) :
    ∃ (Y : Istr P) (φ : X ⟶ Y), IsFrobeniusType (istrPre P) φ ∧
      (istrPre P).degFr φ = n := by
  obtain ⟨B₀, φ₀, hφ, hd⟩ := F.frobDegSurj X.obj n
  exact ⟨⟨B₀, F.isotropicClosed φ₀ X.property⟩, InducedCategory.homMk φ₀,
    (istr_frobType_iff P F _).mpr hφ, hd⟩

include F in
/-- **(v)(b)** の移送。 -/
theorem istr_preStepFactor {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P) φ) :
    ∃ (Z : Istr P) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsCoAngular (istrPre P) β ∧ IsPreStep (istrPre P) β ∧
        IsIsometric (istrPre P) α ∧ IsPreStep (istrPre P) α := by
  obtain ⟨Z₀, β₀, α₀, heq, hβc, hβs, hαi, hαs⟩ := F.preStepFactor φ.hom hφ
  refine ⟨⟨Z₀, F.isotropicClosed β₀ X.property⟩, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, ?_, hβs, hαi, hαs⟩
  exact istr_coAngular P _

include F in
/-- **(v)(c)** の移送。 -/
theorem istr_preStepFactor' {X Y : Istr P} (φ : X ⟶ Y) (hφ : IsPreStep (istrPre P) φ) :
    ∃ (Z : Istr P) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsIsometric (istrPre P) β ∧ IsPreStep (istrPre P) β ∧
        IsCoAngular (istrPre P) α ∧ IsPreStep (istrPre P) α := by
  obtain ⟨Z₀, β₀, α₀, heq, hβi, hβs, hαc, hαs⟩ := F.preStepFactor' φ.hom hφ
  refine ⟨⟨Z₀, F.isotropicClosed β₀ X.property⟩, InducedCategory.homMk β₀,
    InducedCategory.homMk α₀, InducedCategory.hom_ext heq, hβi, hβs, ?_, hαs⟩
  exact istr_coAngular P _

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
theorem istr_isPullBack_to {X Y : Istr P} (g : X ⟶ Y) (h : IsPullBack (istrPre P) g) :
    IsPullBack P g.hom := by
  intro Z
  haveI hZb : IsIso (P.Base (hullMap P F Z)) := (hullMap_spec P F Z).2.1.2
  haveI hZe : Epi (hullMap P F Z) := P.totEpiC _ _ _
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have h1 : (f₁ ≫ g.hom : Z ⟶ Y.obj) = f₂ ≫ g.hom := congrArg Prod.fst hp
    have h2 : P.Base f₁ = P.Base f₂ := congrArg Prod.snd hp
    have e₁ : hullMap P F Z ≫ ((hullHomEquiv P F Z X).symm f₁).hom = f₁ :=
      (hullHomEquiv P F Z X).apply_symm_apply f₁
    have e₂ : hullMap P F Z ≫ ((hullHomEquiv P F Z X).symm f₂).hom = f₂ :=
      (hullHomEquiv P F Z X).apply_symm_apply f₂
    have hgg : (hullHomEquiv P F Z X).symm f₁ = (hullHomEquiv P F Z X).symm f₂ := by
      refine (h (hullIstr P F Z)).1 (Subtype.ext (Prod.ext ?_ ?_))
      · refine InducedCategory.hom_ext ?_
        refine (cancel_epi (hullMap P F Z)).mp ?_
        show hullMap P F Z ≫ (((hullHomEquiv P F Z X).symm f₁).hom ≫ g.hom)
          = hullMap P F Z ≫ (((hullHomEquiv P F Z X).symm f₂).hom ≫ g.hom)
        rw [← Category.assoc, ← Category.assoc, e₁, e₂, h1]
      · show P.Base ((hullHomEquiv P F Z X).symm f₁).hom
          = P.Base ((hullHomEquiv P F Z X).symm f₂).hom
        refine (cancel_epi (P.Base (hullMap P F Z))).mp ?_
        rw [← P.Base_comp, ← P.Base_comp, e₁, e₂, h2]
    rw [← e₁, ← e₂, hgg]
  · rintro ⟨⟨a, b⟩, hab⟩
    have ea : hullMap P F Z ≫ ((hullHomEquiv P F Z Y).symm a).hom = a :=
      (hullHomEquiv P F Z Y).apply_symm_apply a
    have hcond : (istrPre P).Base ((hullHomEquiv P F Z Y).symm a)
        = (inv (P.Base (hullMap P F Z)) ≫ b) ≫ (istrPre P).Base g := by
      show P.Base ((hullHomEquiv P F Z Y).symm a).hom
        = (inv (P.Base (hullMap P F Z)) ≫ b) ≫ P.Base g.hom
      rw [Category.assoc, ← hab, ← ea, P.Base_comp, ← Category.assoc, IsIso.inv_hom_id,
        Category.id_comp]
    obtain ⟨f', hf'⟩ := (h (hullIstr P F Z)).2
      ⟨((hullHomEquiv P F Z Y).symm a, inv (P.Base (hullMap P F Z)) ≫ b), hcond⟩
    have hp := Subtype.ext_iff.mp hf'
    have k1 : (f' ≫ g : hullIstr P F Z ⟶ Y) = (hullHomEquiv P F Z Y).symm a :=
      congrArg Prod.fst hp
    have k2 : P.Base f'.hom = inv (P.Base (hullMap P F Z)) ≫ b := congrArg Prod.snd hp
    refine ⟨hullMap P F Z ≫ f'.hom, Subtype.ext (Prod.ext ?_ ?_)⟩
    · show (hullMap P F Z ≫ f'.hom) ≫ g.hom = a
      rw [Category.assoc, show f'.hom ≫ g.hom = ((hullHomEquiv P F Z Y).symm a).hom from
        congrArg InducedCategory.Hom.hom k1, ea]
    · show P.Base (hullMap P F Z ≫ f'.hom) = b
      rw [P.Base_comp, k2, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

end Istr

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

end ABC3.Found.FrdI
