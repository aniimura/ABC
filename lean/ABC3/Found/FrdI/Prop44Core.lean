import ABC3.Found.FrdI.Prop44Pre

/-!
# [FrdI] Proposition 4.4, (ii) —— `𝒞^birat` の `FrobenioidCore` を埋める

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.85。

原文 (FrdI p.85):
> exercise to check that Cbirat is, in fact, a Frobenioid of group-like type. Moreover, it

★★**要は「`𝒞^birat` の co-angular pre-step はすべて同型」**
(`birat_isIso_of_coaPre_birat`、`Prop44Pre.lean`)である。
★`Definition 1.3` の (iii)(c)・(v)(b)(c)・(vi) はどれも co-angular pre-step についての
条件なので、それが同型に潰れると**ほぼ自明**になる。

## ★本ファイルで埋まる 6 条(測定、2026-08-17)

| 条 | フィールド | 手 |
|---|---|---|
| (v)(b) | `preStepFactor` | ★`φ = 𝟙 ≫ φ`(等長は `𝒞^birat` で自動) |
| (v)(b) | `preStepFactorUniq` | ★`γ := β⁻¹ ≫ β'`(`β`・`β'` は同型) |
| (v)(c) | `preStepFactor'` | ★`φ = φ ≫ 𝟙` |
| (v)(c) | `preStepFactorUniq'` | ★`γ := α ≫ α'⁻¹` |
| (vi) | `faithfulUpToUnits` | ★`α := ψ⁻¹ ≫ φ` |
| (i)(b) | `preStepSpan` | ★`𝒞` のものを押し出す(pre-step は渡る) |

★★**`inv` は使わない** —— `IsIso.out` で逆射を平の射として取り出す
(`lean-rebind-morphisms-clean-types` の規則)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★0. 逆射を平の射として取り出す -/

include P in
/-- ★**co-angular pre-step の逆射**(平の射として)。 -/
theorem birat_coaPre_inv {X Y : BiratCat P G} (f : X ⟶ Y)
    (hc : IsCoAngular (biratPre P G) f) (hs : IsPreStep (biratPre P G) f) :
    ∃ f' : Y ⟶ X, f ≫ f' = 𝟙 X ∧ f' ≫ f = 𝟙 Y :=
  (birat_isIso_of_coaPre_birat f hc hs).out

include P in
/-- ★**同型の Frobenius 次数は 1**(平の射の形で)。 -/
theorem birat_degFr_of_inv {X Y : BiratCat P G} {f : X ⟶ Y} {f' : Y ⟶ X}
    (h : f ≫ f' = 𝟙 X) (hf : (biratPre P G).degFr f = 1) :
    (biratPre P G).degFr f' = 1 := by
  have h1 := congrArg (biratPre P G).degFr h
  rw [(biratPre P G).degFr_comp, hf, mul_one, (biratPre P G).degFr_id] at h1
  exact h1

include P in
/-- ★**同型の底は逆向きの底**(平の射の形で)。 -/
theorem birat_base_of_inv {X Y : BiratCat P G} {f : X ⟶ Y} {f' : Y ⟶ X}
    (h : f' ≫ f = 𝟙 Y) :
    (biratPre P G).Base f' ≫ (biratPre P G).Base f = 𝟙 _ := by
  have h1 := congrArg (biratPre P G).Base h
  rw [(biratPre P G).Base_comp, (biratPre P G).Base_id] at h1
  exact h1

/-! ## ★1. (v)(b)(c) —— pre-step の 2 通りの分解 -/

include P in
theorem birat_preStepFactor {A B : BiratCat P G} (φ : A ⟶ B)
    (hφ : IsPreStep (biratPre P G) φ) :
    ∃ (X : BiratCat P G) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsCoAngular (biratPre P G) β ∧ IsPreStep (biratPre P G) β
        ∧ IsIsometric (biratPre P G) α ∧ IsPreStep (biratPre P G) α :=
  ⟨A, 𝟙 A, φ, (Category.id_comp φ).symm,
    isCoAngular_of_isIso (biratPre P G) (𝟙 A),
    isPreStep_of_isIso (biratPre P G) (𝟙 A), birat_isIsometric φ, hφ⟩

include P in
theorem birat_preStepFactor' {A B : BiratCat P G} (φ : A ⟶ B)
    (hφ : IsPreStep (biratPre P G) φ) :
    ∃ (X : BiratCat P G) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsIsometric (biratPre P G) β ∧ IsPreStep (biratPre P G) β
        ∧ IsCoAngular (biratPre P G) α ∧ IsPreStep (biratPre P G) α :=
  ⟨B, φ, 𝟙 B, (Category.comp_id φ).symm, birat_isIsometric φ, hφ,
    isCoAngular_of_isIso (biratPre P G) (𝟙 B), isPreStep_of_isIso (biratPre P G) (𝟙 B)⟩

/-! ## ★2. (v)(b)(c) の一意性 —— co-angular 部分が同型だから

★★**`γ` は「同型どうしの繋ぎ替え」で作れる。** -/

include P in
theorem birat_preStepFactorUniq {A B : BiratCat P G} (X X' : BiratCat P G)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B) (heq : β ≫ α = β' ≫ α')
    (hβc : IsCoAngular (biratPre P G) β) (hβs : IsPreStep (biratPre P G) β)
    (_hαi : IsIsometric (biratPre P G) α) (_hαs : IsPreStep (biratPre P G) α)
    (hβ'c : IsCoAngular (biratPre P G) β') (hβ's : IsPreStep (biratPre P G) β')
    (_hα'i : IsIsometric (biratPre P G) α') (_hα's : IsPreStep (biratPre P G) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨b, hb1, hb2⟩ := birat_coaPre_inv P G β hβc hβs
  obtain ⟨b', hb'1, hb'2⟩ := birat_coaPre_inv P G β' hβ'c hβ's
  refine ⟨⟨b ≫ β', b' ≫ β, ?_, ?_⟩, ?_, ?_⟩
  · calc (b ≫ β') ≫ (b' ≫ β) = b ≫ ((β' ≫ b') ≫ β) := by simp only [Category.assoc]
      _ = b ≫ (𝟙 A ≫ β) := by rw [hb'1]
      _ = b ≫ β := by rw [Category.id_comp]
      _ = 𝟙 X := hb2
  · calc (b' ≫ β) ≫ (b ≫ β') = b' ≫ ((β ≫ b) ≫ β') := by simp only [Category.assoc]
      _ = b' ≫ (𝟙 A ≫ β') := by rw [hb1]
      _ = b' ≫ β' := by rw [Category.id_comp]
      _ = 𝟙 X' := hb'2
  · calc α' = 𝟙 X' ≫ α' := (Category.id_comp α').symm
      _ = (b' ≫ β') ≫ α' := by rw [hb'2]
      _ = b' ≫ (β' ≫ α') := Category.assoc _ _ _
      _ = b' ≫ (β ≫ α) := by rw [heq]
      _ = (b' ≫ β) ≫ α := (Category.assoc _ _ _).symm
  · calc β' = 𝟙 A ≫ β' := (Category.id_comp β').symm
      _ = (β ≫ b) ≫ β' := by rw [hb1]
      _ = β ≫ (b ≫ β') := Category.assoc _ _ _

include P in
theorem birat_preStepFactorUniq' {A B : BiratCat P G} (X X' : BiratCat P G)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B) (heq : β ≫ α = β' ≫ α')
    (_hβi : IsIsometric (biratPre P G) β) (_hβs : IsPreStep (biratPre P G) β)
    (hαc : IsCoAngular (biratPre P G) α) (hαs : IsPreStep (biratPre P G) α)
    (_hβ'i : IsIsometric (biratPre P G) β') (_hβ's : IsPreStep (biratPre P G) β')
    (hα'c : IsCoAngular (biratPre P G) α') (hα's : IsPreStep (biratPre P G) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  obtain ⟨a, ha1, ha2⟩ := birat_coaPre_inv P G α hαc hαs
  obtain ⟨a', ha'1, ha'2⟩ := birat_coaPre_inv P G α' hα'c hα's
  refine ⟨⟨α ≫ a', α' ≫ a, ?_, ?_⟩, ?_, ?_⟩
  · calc (α ≫ a') ≫ (α' ≫ a) = α ≫ ((a' ≫ α') ≫ a) := by simp only [Category.assoc]
      _ = α ≫ (𝟙 B ≫ a) := by rw [ha'2]
      _ = α ≫ a := by rw [Category.id_comp]
      _ = 𝟙 X := ha1
  · calc (α' ≫ a) ≫ (α ≫ a') = α' ≫ ((a ≫ α) ≫ a') := by simp only [Category.assoc]
      _ = α' ≫ (𝟙 B ≫ a') := by rw [ha2]
      _ = α' ≫ a' := by rw [Category.id_comp]
      _ = 𝟙 X' := ha'1
  · calc α' = α' ≫ 𝟙 B := (Category.comp_id α').symm
      _ = α' ≫ (a ≫ α) := by rw [ha2]
      _ = (α' ≫ a) ≫ α := (Category.assoc _ _ _).symm
  · calc β' = β' ≫ 𝟙 X' := (Category.comp_id β').symm
      _ = β' ≫ (α' ≫ a') := by rw [ha'1]
      _ = (β' ≫ α') ≫ a' := (Category.assoc _ _ _).symm
      _ = (β ≫ α) ≫ a' := by rw [heq]
      _ = β ≫ (α ≫ a') := Category.assoc _ _ _

/-! ## ★3. (vi) —— 単元を除く忠実性

★★co-angular pre-step は同型なので、`α := ψ⁻¹ ≫ φ` が求める単元である。 -/

include P in
theorem birat_faithfulUpToUnits {A B : BiratCat P G} (φ ψ : A ⟶ B)
    (hbe : BaseEquivalent (biratPre P G) φ ψ)
    (_hme : MetricallyEquivalent (biratPre P G) φ ψ)
    (hφc : IsCoAngular (biratPre P G) φ) (hφs : IsPreStep (biratPre P G) φ)
    (hψc : IsCoAngular (biratPre P G) ψ) (hψs : IsPreStep (biratPre P G) ψ) :
    ∃ α : End B, α ∈ OTimes (biratPre P G) B ∧ φ = ψ ≫ (α : B ⟶ B) := by
  obtain ⟨p, hp1, hp2⟩ := birat_coaPre_inv P G φ hφc hφs
  obtain ⟨q, hq1, hq2⟩ := birat_coaPre_inv P G ψ hψc hψs
  refine ⟨q ≫ φ, ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
  · -- ★base-identity
    show (biratPre P G).Base (q ≫ φ) = (biratPre P G).Base (𝟙 B)
    rw [(biratPre P G).Base_comp, hbe, (biratPre P G).Base_id]
    exact birat_base_of_inv P G hq2
  · -- ★linear
    show (biratPre P G).degFr (q ≫ φ) = 1
    rw [(biratPre P G).degFr_comp,
      show (biratPre P G).degFr φ = 1 from hφs.1,
      show (biratPre P G).degFr q = 1 from birat_degFr_of_inv P G hq1 hψs.1, mul_one]
  · -- ★単元
    refine (CategoryTheory.isUnit_iff_isIso ((q ≫ φ : End B) : B ⟶ B)).mpr ?_
    refine ⟨⟨p ≫ ψ, ?_, ?_⟩⟩
    · calc (q ≫ φ) ≫ (p ≫ ψ) = q ≫ ((φ ≫ p) ≫ ψ) := by simp only [Category.assoc]
        _ = q ≫ (𝟙 A ≫ ψ) := by rw [hp1]
        _ = q ≫ ψ := by rw [Category.id_comp]
        _ = 𝟙 B := hq2
    · calc (p ≫ ψ) ≫ (q ≫ φ) = p ≫ ((ψ ≫ q) ≫ φ) := by simp only [Category.assoc]
        _ = p ≫ (𝟙 A ≫ φ) := by rw [hq1]
        _ = p ≫ φ := by rw [Category.id_comp]
        _ = 𝟙 B := hp2
  · -- ★`φ = ψ ≫ (q ≫ φ)`
    calc φ = 𝟙 A ≫ φ := (Category.id_comp φ).symm
      _ = (ψ ≫ q) ≫ φ := by rw [hq1]
      _ = ψ ≫ (q ≫ φ) := Category.assoc _ _ _

/-! ## ★4. ★★★★辞書の残りへ —— 「`𝒞^birat` の射は `[a]⁻¹ ≫ [φ]`」

原文 (FrdI p.84):
> Proposition 1.11, (vii), that there exists a commutative diagram

★★原文は (iv) の co-angular の条を
「**`𝒞^birat` の分解は `𝒞` の分解から来る**([cf. the proof of assertion (i)])」
で片づける。★その「来る」の中身は **`Proposition 1.11, (vii)` を繰り返し当てて
`[·]⁻¹` を左端へ寄せる**ことである。

★本節はその**部品**を 2 つ置く:

| 部品 | 内容 |
|---|---|
| `birat_hom_repr` | `𝒞^birat` の射は `a' ≫ [φ]`(`a'` は co-angular pre-step の逆射) |
| `birat_move_inv` | `[g] ≫ [b]⁻¹ = [g₁]⁻¹ ≫ [k₁]`(`Proposition 1.11, (vii)`) |
-/

include P in
/-- ★★**`𝒞^birat` の射の標準形** —— `f = a' ≫ [φ]`。

★`a'` は co-angular pre-step `a` の**逆射**(平の射として持つ)。 -/
theorem birat_hom_repr {X Y : BiratCat P G} (f : X ⟶ Y) :
    ∃ (A' : C) (a : A' ⟶ biratDown P G X) (φ : A' ⟶ biratDown P G Y)
      (aa : (show BiratCat P G from A') ⟶ X) (a' : X ⟶ (show BiratCat P G from A')),
      IsCoAngular P a ∧ IsPreStep P a ∧ aa = (toBiratCat P G).map a ∧
        aa ≫ a' = 𝟙 _ ∧ a' ≫ aa = 𝟙 X ∧ f = a' ≫ (toBiratCat P G).map φ := by
  obtain ⟨Z, φ, hZφ⟩ := HomBirat.exists_rep f
  obtain ⟨aa, haa⟩ : ∃ aa : (show BiratCat P G from Z.unop.left.obj) ⟶ X,
      aa = (toBiratCat P G).map Z.unop.hom.hom := ⟨_, rfl⟩
  obtain ⟨a', ha1, ha2⟩ : ∃ a' : X ⟶ (show BiratCat P G from Z.unop.left.obj),
      aa ≫ a' = 𝟙 _ ∧ a' ≫ aa = 𝟙 X := by
    rw [haa]
    exact (birat_isIso_of_coaPre Z.unop.hom.hom Z.unop.hom.property.1
      Z.unop.hom.property.2).out
  have hkey : aa ≫ f = (toBiratCat P G).map φ := by
    rw [haa, ← hZφ]
    exact birat_toHom_comp_mk Z.unop.hom.hom Z.unop.hom.property.1
      Z.unop.hom.property.2 φ
  refine ⟨Z.unop.left.obj, Z.unop.hom.hom, φ, aa, a',
    Z.unop.hom.property.1, Z.unop.hom.property.2, haa, ha1, ha2, ?_⟩
  calc f = 𝟙 X ≫ f := (Category.id_comp f).symm
    _ = (a' ≫ aa) ≫ f := by rw [ha2]
    _ = a' ≫ (aa ≫ f) := Category.assoc _ _ _
    _ = a' ≫ (toBiratCat P G).map φ := congrArg (fun t => a' ≫ t) hkey

include P G in
/-- ★★**`[·]⁻¹` を左へ寄せる** —— `Proposition 1.11, (vii)` そのもの。

原文 (FrdI p.84):
> Proposition 1.11, (vii), that there exists a commutative diagram

★`g : A₁ ⟶ U` と co-angular pre-step `b : B₁ ⟶ U` に対し、
`g₁ ≫ g = k₁ ≫ b`(`g₁` は co-angular pre-step)が取れる。
★★これが `𝒞^birat` で `[g] ≫ [b]⁻¹ = [g₁]⁻¹ ≫ [k₁]` を与える。 -/
theorem birat_move_inv {A₁ U B₁ : C} (g : A₁ ⟶ U) (b : B₁ ⟶ U)
    (hbc : IsCoAngular P b) (hbs : IsPreStep P b) :
    ∃ (E : C) (g₁ : E ⟶ A₁) (k₁ : E ⟶ B₁),
      IsCoAngular P g₁ ∧ IsPreStep P g₁ ∧ g₁ ≫ g = k₁ ≫ b := by
  exact prop_1_11_vii P G.core G g b hbc hbs

/-! ## ★5. ★★★★辞書の co-angular の条 —— まず `⟹`(反射)

原文 (FrdI p.85):
> is a co-angular pre-step, so Cbirat

★★**`⟹` は 3 行で落ちる**:
1. `𝒞` の分解を `𝒞^birat` へ押し出す(次数・底・pre-step は渡る、等長は自動)
2. `[φ]` の co-angular 性で **真ん中 `[β]` が同型**
3. ★辞書(`birat_isIso_iff`)で `β` は `𝒞` の **co-angular** pre-step、
   ★★もともと**等長**だったので `Proposition 1.4, (iii)` で**同型**

★★★**要点は 3** —— 「LB-invertible な pre-step は同型」がここで効く。 -/

include P in
/-- ★★★**`[φ]` が co-angular なら `φ` も co-angular**(辞書の `⟹`)。 -/
theorem birat_coAngular_reflect {A B : C} (φ : A ⟶ B)
    (h : IsCoAngular (biratPre P G) ((toBiratCat P G).map φ)) : IsCoAngular P φ := by
  intro X Y γ β α heq hαl hβi hβs hor
  have h1 : (toBiratCat P G).map φ
      = (toBiratCat P G).map γ ≫ (toBiratCat P G).map β ≫ (toBiratCat P G).map α := by
    rw [← (toBiratCat P G).map_comp, ← (toBiratCat P G).map_comp, ← heq]
  have hiso : IsIso ((toBiratCat P G).map β) :=
    h (show BiratCat P G from X) (show BiratCat P G from Y)
      ((toBiratCat P G).map γ) ((toBiratCat P G).map β) ((toBiratCat P G).map α) h1
      ((birat_isLinear_iff α).mpr hαl) (birat_isIsometric _)
      ((birat_isPreStep_iff β).mpr hβs)
      (hor.imp (fun hh => (birat_isBaseIsomorphism_iff α).mpr hh)
        (fun hh => (birat_isBaseIsomorphism_iff γ).mpr hh))
  have hβc : IsCoAngular P β := ((birat_isIso_iff β).mp hiso).1
  exact prop_1_4_iii P G.core β ⟨hβc, hβi⟩ hβs

/-! ## ★6. 同型で挟んでも co-angular は変わらない

★★辞書の `⟸` で、`𝒞^birat` の射を `a' ≫ [φ]` の形に直すときに要る。 -/

/-- ★**同型を前から合成しても co-angular**。 -/
theorem isCoAngular_isoComp {Q : PreFrobenioid C Φ} {A A' B : C} (e : A ⟶ A') (e' : A' ⟶ A)
    (he1 : e ≫ e' = 𝟙 A) (he2 : e' ≫ e = 𝟙 A') {f : A' ⟶ B} (hf : IsCoAngular Q f) :
    IsCoAngular Q (e ≫ f) := by
  intro X Y γ β α heq hαl hβi hβs hor
  refine hf X Y (e' ≫ γ) β α ?_ hαl hβi hβs ?_
  · calc f = 𝟙 A' ≫ f := (Category.id_comp f).symm
      _ = (e' ≫ e) ≫ f := by rw [he2]
      _ = e' ≫ (e ≫ f) := Category.assoc _ _ _
      _ = e' ≫ (γ ≫ β ≫ α) := by rw [heq]
      _ = (e' ≫ γ) ≫ β ≫ α := (Category.assoc _ _ _).symm
  · refine hor.imp id (fun hh => ?_)
    haveI hei : IsIso (Q.Base e) := ⟨⟨Q.Base e', by
        rw [← Q.Base_comp, he1, Q.Base_id], by rw [← Q.Base_comp, he2, Q.Base_id]⟩⟩
    haveI hei' : IsIso (Q.Base e') := ⟨⟨Q.Base e, by
        rw [← Q.Base_comp, he2, Q.Base_id], by rw [← Q.Base_comp, he1, Q.Base_id]⟩⟩
    haveI : IsIso (Q.Base γ) := hh
    show IsIso (Q.Base (e' ≫ γ))
    rw [Q.Base_comp]
    exact IsIso.comp_isIso' hei' hh

/-- ★**同型を後ろから合成しても co-angular**。 -/
theorem isCoAngular_compIso {Q : PreFrobenioid C Φ} {A B B' : C} {f : A ⟶ B}
    (e : B ⟶ B') (e' : B' ⟶ B) (he1 : e ≫ e' = 𝟙 B) (_he2 : e' ≫ e = 𝟙 B')
    (hf : IsCoAngular Q f) : IsCoAngular Q (f ≫ e) := by
  intro X Y γ β α heq hαl hβi hβs hor
  haveI hei : IsIso (Q.Base e) := ⟨⟨Q.Base e', by
      rw [← Q.Base_comp, he1, Q.Base_id], by rw [← Q.Base_comp, _he2, Q.Base_id]⟩⟩
  haveI hei' : IsIso (Q.Base e') := ⟨⟨Q.Base e, by
      rw [← Q.Base_comp, _he2, Q.Base_id], by rw [← Q.Base_comp, he1, Q.Base_id]⟩⟩
  refine hf X Y γ β (α ≫ e') ?_ ?_ hβi hβs ?_
  · calc f = f ≫ 𝟙 B := (Category.comp_id f).symm
      _ = (f ≫ e) ≫ e' := by rw [← he1, Category.assoc]
      _ = (γ ≫ β ≫ α) ≫ e' := by rw [heq]
      _ = γ ≫ β ≫ (α ≫ e') := by simp only [Category.assoc]
  · haveI : IsIso e' := ⟨⟨e, _he2, he1⟩⟩
    show Q.degFr (α ≫ e') = 1
    rw [Q.degFr_comp, show Q.degFr e' = 1 from degFr_of_isIso Q e',
      show Q.degFr α = 1 from hαl, mul_one]
  · refine hor.imp (fun hh => ?_) id
    haveI : IsIso (Q.Base α) := hh
    show IsIso (Q.Base (α ≫ e'))
    rw [Q.Base_comp]
    exact IsIso.comp_isIso' hh hei'

/-! ## ★7. ★★★★★辞書の co-angular の条 —— `⟸`(押し出し)

原文 (FrdI p.85):
> is a co-angular pre-step, so Cbirat

★★**原文の「`𝒞^birat` の分解は `𝒞` の分解から来る」を実体にする。**

| 段 | 内容 |
|---|---|
| 1 | `α = d' ≫ [m]`(標準形)。`β ≫ d'` は pre-step |
| 2 | `β ≫ d' = b' ≫ [k]`。★`k` は `𝒞` の **pre-step** |
| 3 | `γ ≫ b' = c' ≫ [g]` |
| 4 | 合わせて `[c ≫ φ] = [g ≫ k ≫ m]`、忠実性で `c ≫ φ = g ≫ k ≫ m` |
| 5 | ★★`k` を `Definition 1.3, (v), (b)` で `β₀ ≫ α₀`(co-angular / **等長**)に割る |
| 6 | `c ≫ φ` の co-angular 性を `(g ≫ β₀) ≫ α₀ ≫ m` に当てて `α₀` が同型 |
| 7 | よって `k` は co-angular pre-step、`[k]` は同型、`β` も同型 |

★★★**要点は 5** —— 等長でない `k` はそのままでは co-angular 性に食わせられない。 -/

include P in
/-- ★★★★**`φ` が co-angular なら `[φ]` も co-angular**(辞書の `⟸`)。 -/
theorem birat_coAngular_push {A B : C} (φ : A ⟶ B) (hφ : IsCoAngular P φ) :
    IsCoAngular (biratPre P G) ((toBiratCat P G).map φ) := by
  intro U V γ β α heq hαl _hβi hβs hor
  -- ★段 1: `α` の標準形
  obtain ⟨D₁, d, m, dd, d', hdc, hds, hdd, hd1, hd2, hαeq⟩ := birat_hom_repr P G α
  haveI hdiso : IsIso d' := ⟨⟨dd, hd2, hd1⟩⟩
  haveI hddiso : IsIso dd := ⟨⟨d', hd1, hd2⟩⟩
  -- ★段 2: `β ≫ d'` の標準形。★`k` は `𝒞` の pre-step
  have hβd : IsPreStep (biratPre P G) (β ≫ d') :=
    IsPreStep.comp (biratPre P G) hβs (isPreStep_of_isIso (biratPre P G) d')
  obtain ⟨B₂, b, k, bb, b', hbc, hbs, hbb, hb1, hb2, hβeq⟩ := birat_hom_repr P G (β ≫ d')
  haveI hbiso : IsIso b' := ⟨⟨bb, hb2, hb1⟩⟩
  haveI hbbiso : IsIso bb := ⟨⟨b', hb1, hb2⟩⟩
  have hkps : IsPreStep P k := by
    refine (birat_isPreStep_iff (P := P) (G := G) k).mp ?_
    have h1 : IsPreStep (biratPre P G) (b' ≫ (toBiratCat P G).map k) := hβeq ▸ hβd
    have h3 : bb ≫ (b' ≫ (toBiratCat P G).map k) = (toBiratCat P G).map k :=
      ((Category.assoc _ _ _).symm.trans
        (congrArg (fun t => t ≫ (toBiratCat P G).map k) hb1)).trans (Category.id_comp _)
    exact h3 ▸ IsPreStep.comp (biratPre P G)
      (isPreStep_of_isIso (biratPre P G) bb) h1
  -- ★段 3: `γ ≫ b'` の標準形
  obtain ⟨C₂, c, g, cc, c', hcc, hcs, hccdef, hc1, hc2, hγeq⟩ := birat_hom_repr P G (γ ≫ b')
  haveI hc'iso : IsIso c' := ⟨⟨cc, hc2, hc1⟩⟩
  haveI hcciso : IsIso cc := ⟨⟨c', hc1, hc2⟩⟩
  -- ★段 4: `𝒞` の等式
  have hkey0 : (toBiratCat P G).map φ
      = c' ≫ ((toBiratCat P G).map g
        ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m)) := by
    calc (toBiratCat P G).map φ = γ ≫ (β ≫ α) := heq
      _ = γ ≫ ((β ≫ d') ≫ (toBiratCat P G).map m) := by
          rw [hαeq]; simp only [Category.assoc]; try rfl
      _ = γ ≫ ((b' ≫ (toBiratCat P G).map k) ≫ (toBiratCat P G).map m) := by rw [hβeq]; try rfl
      _ = (γ ≫ b') ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m) := by
          simp only [Category.assoc]; try rfl
      _ = (c' ≫ (toBiratCat P G).map g)
            ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m) := by rw [hγeq]; try rfl
      _ = c' ≫ ((toBiratCat P G).map g
            ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m)) := Category.assoc _ _ _
  haveI := toBiratCat_faithful (P := P) (G := G)
  have hkey : c ≫ φ = g ≫ (k ≫ m) := by
    refine (toBiratCat P G).map_injective ?_
    have e1 : (toBiratCat P G).map (c ≫ φ) = cc ≫ (toBiratCat P G).map φ := by
      rw [(toBiratCat P G).map_comp, ← hccdef]; try rfl
    have e2 : cc ≫ (toBiratCat P G).map φ
        = cc ≫ (c' ≫ ((toBiratCat P G).map g
          ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m))) :=
      congrArg (fun t => cc ≫ t) hkey0
    have e3 : cc ≫ (c' ≫ ((toBiratCat P G).map g
          ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m)))
        = (toBiratCat P G).map g ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m) :=
      ((Category.assoc _ _ _).symm.trans
        (congrArg (fun t => t ≫ ((toBiratCat P G).map g
          ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m))) hc1)).trans
        (Category.id_comp _)
    have e4 : (toBiratCat P G).map g ≫ ((toBiratCat P G).map k ≫ (toBiratCat P G).map m)
        = (toBiratCat P G).map (g ≫ (k ≫ m)) :=
      (congrArg (fun t => (toBiratCat P G).map g ≫ t)
        ((toBiratCat P G).map_comp k m).symm).trans
        ((toBiratCat P G).map_comp g (k ≫ m)).symm
    exact ((e1.trans e2).trans e3).trans e4
  -- ★段 5: `k` を等長部分と co-angular 部分に割る
  obtain ⟨K₀, β₀, α₀, hkeq, hβ₀c, hβ₀s, hα₀i, hα₀s⟩ := G.core.preStepFactor k hkps
  -- ★`m` は linear
  have hml : IsLinear P m := by
    have h2 : (biratPre P G).degFr ((toBiratCat P G).map m) * (biratPre P G).degFr d' = 1 :=
      ((biratPre P G).degFr_comp d' ((toBiratCat P G).map m)).symm.trans
        (show (biratPre P G).degFr (d' ≫ (toBiratCat P G).map m) = 1 from hαeq ▸ hαl)
    rw [show (biratPre P G).degFr d' = 1 from degFr_of_isIso (biratPre P G) d',
      mul_one] at h2
    exact (birat_isLinear_iff (P := P) (G := G) m).mp h2
  -- ★段 6: `c ≫ φ` の co-angular 性
  have hcφ : IsCoAngular P (c ≫ φ) := IsCoAngular.comp P G.core hcc hφ
  have heq6 : c ≫ φ = (g ≫ β₀) ≫ α₀ ≫ m := by
    rw [hkey, hkeq]
    calc g ≫ ((β₀ ≫ α₀) ≫ m) = g ≫ (β₀ ≫ (α₀ ≫ m)) :=
          congrArg (fun t => g ≫ t) (Category.assoc _ _ _)
      _ = (g ≫ β₀) ≫ (α₀ ≫ m) := (Category.assoc _ _ _).symm
  have hor6 : IsBaseIsomorphism P m ∨ IsBaseIsomorphism P (g ≫ β₀) := by
    refine hor.imp (fun hh => ?_) (fun hh => ?_)
    · -- `Base α` 同型 ⟹ `Base m` 同型
      haveI : IsIso ((biratPre P G).Base α) := hh
      haveI hbd : IsIso ((biratPre P G).Base d') := isBaseIsomorphism_of_isIso _ d'
      refine (birat_isBaseIsomorphism_iff (P := P) (G := G) m).mp ?_
      haveI h1 : IsIso ((biratPre P G).Base d' ≫ (biratPre P G).Base
          ((toBiratCat P G).map m)) := by
        rw [← (biratPre P G).Base_comp, ← hαeq]; exact hh
      exact IsIso.of_isIso_comp_left ((biratPre P G).Base d')
        ((biratPre P G).Base ((toBiratCat P G).map m))
    · -- `Base γ` 同型 ⟹ `Base (g ≫ β₀)` 同型
      haveI : IsIso ((biratPre P G).Base γ) := hh
      haveI hgb : IsIso ((biratPre P G).Base b') := isBaseIsomorphism_of_isIso _ b'
      haveI hgc : IsIso ((biratPre P G).Base c') := isBaseIsomorphism_of_isIso _ c'
      have e1 : (biratPre P G).Base c' ≫ (biratPre P G).Base ((toBiratCat P G).map g)
          = (biratPre P G).Base γ ≫ (biratPre P G).Base b' := by
        rw [← (biratPre P G).Base_comp, ← (biratPre P G).Base_comp, hγeq]; try rfl
      haveI h1 : IsIso ((biratPre P G).Base c' ≫ (biratPre P G).Base
          ((toBiratCat P G).map g)) := by
        rw [e1]; exact IsIso.comp_isIso' hh hgb
      haveI hgg : IsIso ((biratPre P G).Base ((toBiratCat P G).map g)) :=
        IsIso.of_isIso_comp_left ((biratPre P G).Base c')
          ((biratPre P G).Base ((toBiratCat P G).map g))
      haveI : IsIso (P.Base g) := (birat_isBaseIsomorphism_iff (P := P) (G := G) g).mp hgg
      show IsIso (P.Base (g ≫ β₀))
      rw [P.Base_comp]
      exact IsIso.comp_isIso' inferInstance hβ₀s.2
  haveI hα₀iso : IsIso α₀ := hcφ K₀ _ (g ≫ β₀) α₀ m heq6 hml hα₀i hα₀s hor6
  -- ★段 7: `k` は co-angular pre-step、よって `[k]` は同型
  haveI hkiso : IsIso ((toBiratCat P G).map k) := by
    refine birat_isIso_of_coaPre (P := P) (G := G) k ?_ hkps
    rw [hkeq]
    exact IsCoAngular.comp P G.core hβ₀c (isCoAngular_of_isIso P α₀)
  -- ★`β = (b' ≫ [k]) ≫ dd`
  have hβfin : β = (b' ≫ (toBiratCat P G).map k) ≫ dd := by
    calc β = β ≫ 𝟙 V := (Category.comp_id β).symm
      _ = β ≫ (d' ≫ dd) := by rw [hd2]; try rfl
      _ = (β ≫ d') ≫ dd := (Category.assoc _ _ _).symm
      _ = (b' ≫ (toBiratCat P G).map k) ≫ dd := by rw [hβeq]; try rfl
  rw [hβfin]
  exact IsIso.comp_isIso' (IsIso.comp_isIso' hbiso hkiso) hddiso

/-- ★★★★**[FrdI] Proposition 4.4, (iv) の co-angular の条** —— 辞書。 -/
theorem birat_isCoAngular_iff {A B : C} (φ : A ⟶ B) :
    IsCoAngular (biratPre P G) ((toBiratCat P G).map φ) ↔ IsCoAngular P φ :=
  ⟨birat_coAngular_reflect P G φ, birat_coAngular_push P G φ⟩

end ABC3.Found.FrdI
