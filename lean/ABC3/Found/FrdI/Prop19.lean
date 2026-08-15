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
