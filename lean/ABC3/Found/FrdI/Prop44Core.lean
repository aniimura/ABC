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

/-! ## ★8. 辞書の残り —— Frobenius 型、および「代表元で判定できる」

原文 (FrdI p.83):
> morphism of a given Frobenius degree; isometry; pre-step; base-isomorphism) of

★★`𝒞^birat` では**等長が自動**なので、
`Frobenius 型 = co-angular ∧ 底同型` に潰れる ——
これが原文の「Frobenius 型 ⟺ **co-angular base-isomorphism**」である。 -/

include P in
/-- ★★★**[FrdI] Proposition 4.4, (iv) の「Frobenius 型」の条**。 -/
theorem birat_isFrobeniusType_iff {A B : C} (φ : A ⟶ B) :
    IsFrobeniusType (biratPre P G) ((toBiratCat P G).map φ)
      ↔ (IsCoAngular P φ ∧ IsBaseIsomorphism P φ) := by
  constructor
  · rintro ⟨⟨hc, _⟩, hb⟩
    exact ⟨(birat_isCoAngular_iff P G φ).mp hc,
      (birat_isBaseIsomorphism_iff (P := P) (G := G) φ).mp hb⟩
  · rintro ⟨hc, hb⟩
    exact ⟨⟨(birat_isCoAngular_iff P G φ).mpr hc, birat_isIsometric _⟩,
      (birat_isBaseIsomorphism_iff (P := P) (G := G) φ).mpr hb⟩

include P in
/-- ★★**co-angular 性は代表元で判定できる** —— `f = a' ≫ [φ]` のとき
`f` が co-angular ⟺ `φ` が `𝒞` で co-angular。 -/
theorem birat_isCoAngular_repr {X Y : BiratCat P G} (f : X ⟶ Y)
    {A' : C} (φ : A' ⟶ biratDown P G Y)
    (aa : (show BiratCat P G from A') ⟶ X) (a' : X ⟶ (show BiratCat P G from A'))
    (h1 : aa ≫ a' = 𝟙 _) (h2 : a' ≫ aa = 𝟙 X)
    (hf : f = a' ≫ (toBiratCat P G).map φ) :
    IsCoAngular (biratPre P G) f ↔ IsCoAngular P φ := by
  have hmap : (toBiratCat P G).map φ = aa ≫ f := by
    rw [hf]
    exact (((Category.assoc _ _ _).symm.trans
      (congrArg (fun t => t ≫ (toBiratCat P G).map φ) h1)).trans
      (Category.id_comp _)).symm
  constructor
  · intro hc
    refine (birat_isCoAngular_iff P G φ).mp ?_
    rw [hmap]
    exact isCoAngular_isoComp aa a' h1 h2 hc
  · intro hc
    rw [hf]
    exact isCoAngular_isoComp a' aa h2 h1 ((birat_isCoAngular_iff P G φ).mpr hc)

/-! ## ★9. (iii)(a) —— co-angular は合成で閉じる

原文 (FrdI p.24):
> (iii) (Surjectivity to the Divisor Monoid via Co-angular Morphisms) (a) The

★★**辞書を 2 度使う**:
1. `ψ₀`・`φ₀` は `𝒞` で co-angular(辞書)
2. `Proposition 1.11, (vii)` で `g₁ ≫ ψ₀ = k₁ ≫ b`。
   ★`g₁ ≫ ψ₀` が co-angular なので `k₁ ≫ b` も、
   ★★`[b]` は同型なので **`k₁` も co-angular**(辞書を往復)
3. `ψ ≫ φ = (a' ≫ g') ≫ [k₁ ≫ φ₀]` で、`k₁ ≫ φ₀` は co-angular -/

include P in
/-- ★★★**[FrdI] Definition 1.3, (iii), (a)** for `𝒞^birat`。 -/
theorem birat_coAngularComp {A B E : BiratCat P G} (ψ : A ⟶ B) (φ : B ⟶ E)
    (hψ : IsCoAngular (biratPre P G) ψ) (hφ : IsCoAngular (biratPre P G) φ) :
    IsCoAngular (biratPre P G) (ψ ≫ φ) := by
  obtain ⟨A₁, a, ψ₀, aa, a', hac, has, haa, ha1, ha2, hψeq⟩ := birat_hom_repr P G ψ
  obtain ⟨B₁, b, φ₀, bb, b', hbc, hbs, hbb, hb1, hb2, hφeq⟩ := birat_hom_repr P G φ
  have hψ₀c : IsCoAngular P ψ₀ :=
    (birat_isCoAngular_repr P G ψ ψ₀ aa a' ha1 ha2 hψeq).mp hψ
  have hφ₀c : IsCoAngular P φ₀ :=
    (birat_isCoAngular_repr P G φ φ₀ bb b' hb1 hb2 hφeq).mp hφ
  -- ★`Proposition 1.11, (vii)`
  obtain ⟨F₁, g₁, k₁, hg₁c, hg₁s, hsq⟩ := birat_move_inv P G ψ₀ b hbc hbs
  obtain ⟨gg, hgg⟩ : ∃ gg : (show BiratCat P G from F₁) ⟶ (show BiratCat P G from A₁),
      gg = (toBiratCat P G).map g₁ := ⟨_, rfl⟩
  obtain ⟨g', hg'1, hg'2⟩ : ∃ g' : (show BiratCat P G from A₁) ⟶ (show BiratCat P G from F₁),
      gg ≫ g' = 𝟙 _ ∧ g' ≫ gg = 𝟙 _ := by
    rw [hgg]
    exact (birat_isIso_of_coaPre (P := P) (G := G) g₁ hg₁c hg₁s).out
  -- ★`k₁` は co-angular
  have hk₁c : IsCoAngular P k₁ := by
    have h1 : IsCoAngular P (k₁ ≫ b) := by
      rw [← hsq]
      exact IsCoAngular.comp P G.core hg₁c hψ₀c
    have h2 : IsCoAngular (biratPre P G) ((toBiratCat P G).map (k₁ ≫ b)) :=
      (birat_isCoAngular_iff P G _).mpr h1
    have h3 : (toBiratCat P G).map (k₁ ≫ b) ≫ b' = (toBiratCat P G).map k₁ := by
      rw [(toBiratCat P G).map_comp, ← hbb]
      exact ((Category.assoc _ _ _).trans
        (congrArg (fun t => (toBiratCat P G).map k₁ ≫ t) hb1)).trans (Category.comp_id _)
    refine (birat_isCoAngular_iff P G k₁).mp ?_
    rw [← h3]
    exact isCoAngular_compIso b' bb hb2 hb1 h2
  -- ★`[ψ₀] ≫ b' = g' ≫ [k₁]`
  have e0 : gg ≫ (toBiratCat P G).map ψ₀ = (toBiratCat P G).map k₁ ≫ bb := by
    rw [hgg, hbb]
    exact ((toBiratCat P G).map_comp g₁ ψ₀).symm.trans
      ((congrArg (toBiratCat P G).map hsq).trans ((toBiratCat P G).map_comp k₁ b))
  have hmov : (toBiratCat P G).map ψ₀ ≫ b' = g' ≫ (toBiratCat P G).map k₁ := by
    have e2 : ((toBiratCat P G).map k₁ ≫ bb) ≫ b' = (toBiratCat P G).map k₁ :=
      ((Category.assoc _ _ _).trans
        (congrArg (fun t => (toBiratCat P G).map k₁ ≫ t) hb1)).trans (Category.comp_id _)
    have s1 : g' ≫ (toBiratCat P G).map k₁
        = g' ≫ (((toBiratCat P G).map k₁ ≫ bb) ≫ b') :=
      congrArg (fun t => g' ≫ t) e2.symm
    have s2 : g' ≫ (((toBiratCat P G).map k₁ ≫ bb) ≫ b')
        = g' ≫ ((gg ≫ (toBiratCat P G).map ψ₀) ≫ b') :=
      congrArg (fun t => g' ≫ (t ≫ b')) e0.symm
    have s3 : g' ≫ ((gg ≫ (toBiratCat P G).map ψ₀) ≫ b')
        = (g' ≫ gg) ≫ ((toBiratCat P G).map ψ₀ ≫ b') := by
      simp only [Category.assoc]; try rfl
    have s4 : (g' ≫ gg) ≫ ((toBiratCat P G).map ψ₀ ≫ b')
        = (toBiratCat P G).map ψ₀ ≫ b' :=
      (congrArg (fun t => t ≫ ((toBiratCat P G).map ψ₀ ≫ b')) hg'2).trans
        (Category.id_comp _)
    exact ((s1.trans s2).trans (s3.trans s4)).symm
  -- ★`ψ ≫ φ = (a' ≫ g') ≫ [k₁ ≫ φ₀]`
  have hfin : ψ ≫ φ = (a' ≫ g') ≫ (toBiratCat P G).map (k₁ ≫ φ₀) := by
    have e3 : (toBiratCat P G).map (k₁ ≫ φ₀)
        = (toBiratCat P G).map k₁ ≫ (toBiratCat P G).map φ₀ :=
      (toBiratCat P G).map_comp k₁ φ₀
    rw [hψeq, hφeq, e3]
    have t1 : (a' ≫ (toBiratCat P G).map ψ₀) ≫ (b' ≫ (toBiratCat P G).map φ₀)
        = a' ≫ (((toBiratCat P G).map ψ₀ ≫ b') ≫ (toBiratCat P G).map φ₀) := by
      simp only [Category.assoc]; try rfl
    have t2 : a' ≫ (((toBiratCat P G).map ψ₀ ≫ b') ≫ (toBiratCat P G).map φ₀)
        = a' ≫ ((g' ≫ (toBiratCat P G).map k₁) ≫ (toBiratCat P G).map φ₀) :=
      congrArg (fun t => a' ≫ (t ≫ (toBiratCat P G).map φ₀)) hmov
    have t3 : a' ≫ ((g' ≫ (toBiratCat P G).map k₁) ≫ (toBiratCat P G).map φ₀)
        = (a' ≫ g') ≫ ((toBiratCat P G).map k₁ ≫ (toBiratCat P G).map φ₀) := by
      simp only [Category.assoc]; try rfl
    exact (t1.trans t2).trans t3
  have hi1 : (gg ≫ aa) ≫ (a' ≫ g') = 𝟙 _ := by
    have u1 : (gg ≫ aa) ≫ (a' ≫ g') = gg ≫ ((aa ≫ a') ≫ g') := by
      simp only [Category.assoc]; try rfl
    have u2 : gg ≫ ((aa ≫ a') ≫ g') = gg ≫ (𝟙 _ ≫ g') :=
      congrArg (fun t => gg ≫ (t ≫ g')) ha1
    have u3 : gg ≫ (𝟙 _ ≫ g') = gg ≫ g' :=
      congrArg (fun t => gg ≫ t) (Category.id_comp g')
    exact ((u1.trans u2).trans u3).trans hg'1
  have hi2 : (a' ≫ g') ≫ (gg ≫ aa) = 𝟙 A := by
    have u1 : (a' ≫ g') ≫ (gg ≫ aa) = a' ≫ ((g' ≫ gg) ≫ aa) := by
      simp only [Category.assoc]; try rfl
    have u2 : a' ≫ ((g' ≫ gg) ≫ aa) = a' ≫ (𝟙 _ ≫ aa) :=
      congrArg (fun t => a' ≫ (t ≫ aa)) hg'2
    have u3 : a' ≫ (𝟙 _ ≫ aa) = a' ≫ aa :=
      congrArg (fun t => a' ≫ t) (Category.id_comp aa)
    exact ((u1.trans u2).trans u3).trans ha2
  refine (birat_isCoAngular_repr P G (ψ ≫ φ) (k₁ ≫ φ₀) (gg ≫ aa) (a' ≫ g')
    hi1 hi2 hfin).mpr ?_
  exact IsCoAngular.comp P G.core hk₁c hφ₀c

/-! ## ★10. `𝒞` から押し出すだけの 3 条 —— (i)(a)(b)・(ii)

★★底(`𝒟`)も次数も `𝒞^birat` で変わらないので、`𝒞` の証人をそのまま送る。
★Frobenius 型が渡ることは辞書(`birat_isFrobeniusType_iff`)による。 -/

include P in
/-- ★**(i)(b)** 底の同型は pre-step の span で持ち上がる。 -/
theorem birat_preStepSpan (A B : BiratCat P G)
    (α : ((biratPre P G).toElem.obj A).base ⟶ ((biratPre P G).toElem.obj B).base)
    (hα : IsIso α) :
    ∃ (X : BiratCat P G) (φ : X ⟶ A) (ψ : X ⟶ B) (hφ : IsPreStep (biratPre P G) φ),
      IsPreStep (biratPre P G) ψ ∧
        α = @inv _ _ _ _ ((biratPre P G).Base φ) hφ.2 ≫ (biratPre P G).Base ψ := by
  obtain ⟨X, φ, ψ, hφ, hψ, heq⟩ := G.core.preStepSpan (biratDown P G A) (biratDown P G B) α hα
  refine ⟨show BiratCat P G from X, (toBiratCat P G).map φ, (toBiratCat P G).map ψ,
    (birat_isPreStep_iff (P := P) (G := G) φ).mpr hφ,
    (birat_isPreStep_iff (P := P) (G := G) ψ).mpr hψ, ?_⟩
  haveI : IsIso ((biratPre P G).Base ((toBiratCat P G).map φ)) :=
    ((birat_isPreStep_iff (P := P) (G := G) φ).mpr hφ).2
  refine (IsIso.eq_inv_comp ((biratPre P G).Base ((toBiratCat P G).map φ))).mpr ?_
  have e1 : (biratPre P G).Base ((toBiratCat P G).map φ) = P.Base φ :=
    biratBase_toHomBirat φ
  have e2 : (biratPre P G).Base ((toBiratCat P G).map ψ) = P.Base ψ :=
    biratBase_toHomBirat ψ
  rw [e1]
  haveI : IsIso (P.Base φ) := hφ.2
  exact ((IsIso.eq_inv_comp (P.Base φ)).mp heq).trans e2.symm

include P in
/-- ★**(ii)** 各次数の Frobenius 型射が存在する。 -/
theorem birat_frobDegSurj (A : BiratCat P G) (n : ℕ+) :
    ∃ (B : BiratCat P G) (φ : A ⟶ B),
      IsFrobeniusType (biratPre P G) φ ∧ (biratPre P G).degFr φ = n := by
  obtain ⟨B, φ, hft, hdeg⟩ := G.core.frobDegSurj (biratDown P G A) n
  exact ⟨show BiratCat P G from B, (toBiratCat P G).map φ,
    (birat_isFrobeniusType_iff P G φ).mpr ⟨hft.1.1, hft.2⟩,
    (birat_degFr_iff (P := P) (G := G) φ n).mpr hdeg⟩

include P in
/-- ★**Frobenius-trivial 性は `𝒞^birat` へ渡る** —— `ζ` を関手で送るだけ。 -/
theorem birat_isFrobeniusTrivial (A : C) (h : IsFrobeniusTrivial P A) :
    IsFrobeniusTrivial (biratPre P G) (show BiratCat P G from A) := by
  obtain ⟨ζ, hζd, hζb⟩ := h
  refine ⟨{ toFun := fun n => (toBiratCat P G).map ((ζ n : A ⟶ A))
            map_one' := by
              show (toBiratCat P G).map ((ζ 1 : A ⟶ A)) = 𝟙 _
              rw [show (ζ 1 : A ⟶ A) = 𝟙 A from
                congrArg (fun t : End A => (t : A ⟶ A)) ζ.map_one]
              exact (toBiratCat P G).map_id A
            map_mul' := fun x y => by
              show (toBiratCat P G).map ((ζ (x * y) : A ⟶ A))
                = (toBiratCat P G).map ((ζ y : A ⟶ A))
                  ≫ (toBiratCat P G).map ((ζ x : A ⟶ A))
              rw [show (ζ (x * y) : A ⟶ A) = (ζ y : A ⟶ A) ≫ (ζ x : A ⟶ A) from
                congrArg (fun t : End A => (t : A ⟶ A)) (ζ.map_mul x y)]
              exact (toBiratCat P G).map_comp _ _ }, ?_, ?_⟩
  · intro n
    exact (birat_degFr_iff (P := P) (G := G) ((ζ n : A ⟶ A)) n).mpr (hζd n)
  · intro n
    refine ⟨?_, ?_⟩
    · show (biratPre P G).Base ((toBiratCat P G).map ((ζ n : A ⟶ A)))
        = (biratPre P G).Base (𝟙 _)
      have h1 : (biratPre P G).Base ((toBiratCat P G).map ((ζ n : A ⟶ A)))
          = P.Base ((ζ n : A ⟶ A)) := biratBase_toHomBirat ((ζ n : A ⟶ A))
      rw [h1, (biratPre P G).Base_id]
      exact ((hζb n).1).trans (P.Base_id A)
    · exact (birat_isFrobeniusType_iff P G ((ζ n : A ⟶ A))).mpr
        ⟨(hζb n).2.1.1, (hζb n).2.2⟩

include P in
/-- ★**(i)(a)** 底への全射性。 -/
theorem birat_baseSurj (Y : D) :
    ∃ A : BiratCat P G, IsFrobeniusTrivial (biratPre P G) A ∧
      Nonempty (((biratPre P G).toElem.obj A).base ≅ Y) := by
  obtain ⟨A, hft, he⟩ := G.core.baseSurj Y
  exact ⟨show BiratCat P G from A, birat_isFrobeniusTrivial P G A hft, he⟩

/-! ## ★11. (iii)(c) の全単射 —— co-angular pre-step が同型だから共役で済む

原文 (FrdI p.24):
> depends only [among the bijections induced by the various co-angular pre-steps

★★`𝒞^birat` では co-angular pre-step が**同型**なので、
`𝒪^▷` の対応は**共役** `β = φ⁻¹ ≫ α ≫ φ` そのものである。 -/

include P in
/-- ★★**(iii)(c) の順方向**。 -/
theorem birat_otriFwd {A B : BiratCat P G} (φ : A ⟶ B)
    (hc : IsCoAngular (biratPre P G) φ) (hs : IsPreStep (biratPre P G) φ)
    (α : End A) (hα : α ∈ OTri (biratPre P G) A) :
    ∃! β : End B, β ∈ OTri (biratPre P G) B ∧
      ((φ ≫ (β : B ⟶ B) : A ⟶ B) = ((α : A ⟶ A) ≫ φ)) := by
  obtain ⟨p, hp1, hp2⟩ := birat_coaPre_inv P G φ hc hs
  have hpd : (biratPre P G).degFr p = 1 := birat_degFr_of_inv P G hp1 hs.1
  have hpb : (biratPre P G).Base p ≫ (biratPre P G).Base φ = 𝟙 _ :=
    birat_base_of_inv P G hp2
  refine ⟨p ≫ ((α : A ⟶ A) ≫ φ), ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
  · -- ★base-identity
    show (biratPre P G).Base (p ≫ ((α : A ⟶ A) ≫ φ)) = (biratPre P G).Base (𝟙 B)
    rw [(biratPre P G).Base_comp, (biratPre P G).Base_comp,
      show (biratPre P G).Base (α : A ⟶ A) = (biratPre P G).Base (𝟙 A) from hα.1,
      (biratPre P G).Base_id, Category.id_comp, (biratPre P G).Base_id]
    exact hpb
  · -- ★linear
    show (biratPre P G).degFr (p ≫ ((α : A ⟶ A) ≫ φ)) = 1
    rw [(biratPre P G).degFr_comp, (biratPre P G).degFr_comp,
      show (biratPre P G).degFr φ = 1 from hs.1,
      show (biratPre P G).degFr (α : A ⟶ A) = 1 from hα.2, hpd, one_mul, one_mul]
  · -- ★等式
    have e1 : φ ≫ (p ≫ ((α : A ⟶ A) ≫ φ)) = (φ ≫ p) ≫ ((α : A ⟶ A) ≫ φ) :=
      (Category.assoc _ _ _).symm
    exact e1.trans ((congrArg (fun t => t ≫ ((α : A ⟶ A) ≫ φ)) hp1).trans
      (Category.id_comp _))
  · -- ★一意性
    rintro β' ⟨-, hβ'⟩
    have e2 : (β' : B ⟶ B) = 𝟙 B ≫ (β' : B ⟶ B) := (Category.id_comp _).symm
    have e3 : 𝟙 B ≫ (β' : B ⟶ B) = (p ≫ φ) ≫ (β' : B ⟶ B) :=
      congrArg (fun t => t ≫ (β' : B ⟶ B)) hp2.symm
    have e4 : (p ≫ φ) ≫ (β' : B ⟶ B) = p ≫ (φ ≫ (β' : B ⟶ B)) := Category.assoc _ _ _
    have e5 : p ≫ (φ ≫ (β' : B ⟶ B)) = p ≫ ((α : A ⟶ A) ≫ φ) :=
      congrArg (fun t => p ≫ t) hβ'
    exact ((e2.trans e3).trans (e4.trans e5))

include P in
/-- ★★**(iii)(c) の逆方向**。 -/
theorem birat_otriBwd {A B : BiratCat P G} (φ : A ⟶ B)
    (hc : IsCoAngular (biratPre P G) φ) (hs : IsPreStep (biratPre P G) φ)
    (β : End B) (hβ : β ∈ OTri (biratPre P G) B) :
    ∃! α : End A, α ∈ OTri (biratPre P G) A ∧
      ((φ ≫ (β : B ⟶ B) : A ⟶ B) = ((α : A ⟶ A) ≫ φ)) := by
  obtain ⟨p, hp1, hp2⟩ := birat_coaPre_inv P G φ hc hs
  have hpd : (biratPre P G).degFr p = 1 := birat_degFr_of_inv P G hp1 hs.1
  have hpb : (biratPre P G).Base φ ≫ (biratPre P G).Base p = 𝟙 _ :=
    birat_base_of_inv P G hp1
  refine ⟨φ ≫ ((β : B ⟶ B) ≫ p), ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
  · show (biratPre P G).Base (φ ≫ ((β : B ⟶ B) ≫ p)) = (biratPre P G).Base (𝟙 A)
    rw [(biratPre P G).Base_comp, (biratPre P G).Base_comp,
      show (biratPre P G).Base (β : B ⟶ B) = (biratPre P G).Base (𝟙 B) from hβ.1,
      (biratPre P G).Base_id, Category.id_comp, (biratPre P G).Base_id]
    exact hpb
  · show (biratPre P G).degFr (φ ≫ ((β : B ⟶ B) ≫ p)) = 1
    rw [(biratPre P G).degFr_comp, (biratPre P G).degFr_comp,
      show (biratPre P G).degFr φ = 1 from hs.1,
      show (biratPre P G).degFr (β : B ⟶ B) = 1 from hβ.2, hpd, one_mul, mul_one]
  · -- ★等式: `φ ≫ β = (φ ≫ β ≫ p) ≫ φ`
    have e1 : (φ ≫ ((β : B ⟶ B) ≫ p)) ≫ φ = φ ≫ (((β : B ⟶ B) ≫ p) ≫ φ) :=
      Category.assoc _ _ _
    have e2 : ((β : B ⟶ B) ≫ p) ≫ φ = (β : B ⟶ B) ≫ (p ≫ φ) := Category.assoc _ _ _
    have e3 : (β : B ⟶ B) ≫ (p ≫ φ) = (β : B ⟶ B) ≫ 𝟙 B :=
      congrArg (fun t => (β : B ⟶ B) ≫ t) hp2
    exact (e1.trans ((congrArg (fun t => φ ≫ t) (e2.trans e3)).trans
      (congrArg (fun t => φ ≫ t) (Category.comp_id _)))).symm
  · rintro α' ⟨-, hα'⟩
    have e2 : (α' : A ⟶ A) = (α' : A ⟶ A) ≫ 𝟙 A := (Category.comp_id _).symm
    have e3 : (α' : A ⟶ A) ≫ 𝟙 A = (α' : A ⟶ A) ≫ (φ ≫ p) :=
      congrArg (fun t => (α' : A ⟶ A) ≫ t) hp1.symm
    have e4 : (α' : A ⟶ A) ≫ (φ ≫ p) = ((α' : A ⟶ A) ≫ φ) ≫ p :=
      (Category.assoc _ _ _).symm
    have e5 : ((α' : A ⟶ A) ≫ φ) ≫ p = (φ ≫ (β : B ⟶ B)) ≫ p :=
      congrArg (fun t => t ≫ p) hα'.symm
    exact ((e2.trans e3).trans (e4.trans e5)).trans (Category.assoc _ _ _)

/-! ## ★12. (v)(a) —— pre-step は mono

★★代表元に落として `𝒞` の「pre-step は mono」を使う。 -/

include P in
theorem birat_preStepMono {A B : BiratCat P G} (f : A ⟶ B)
    (hf : IsPreStep (biratPre P G) f) : Mono f := by
  obtain ⟨A₁, a, φ, aa, a', hac, has, haa, ha1, ha2, hfeq⟩ := birat_hom_repr P G f
  have hφs : IsPreStep P φ := by
    refine (birat_isPreStep_iff (P := P) (G := G) φ).mp ?_
    have e1 : (toBiratCat P G).map φ = aa ≫ f := by
      rw [hfeq]
      exact (((Category.assoc _ _ _).symm.trans
        (congrArg (fun t => t ≫ (toBiratCat P G).map φ) ha1)).trans
        (Category.id_comp _)).symm
    rw [e1]
    exact IsPreStep.comp (biratPre P G)
      (haa ▸ (birat_isPreStep_iff (P := P) (G := G) a).mpr has) hf
  haveI := toBiratCat_faithful (P := P) (G := G)
  haveI hmono : Mono φ := G.core.preStepMono φ hφs
  haveI hiso : IsIso a' := ⟨⟨aa, ha2, ha1⟩⟩
  haveI hmono2 : Mono ((toBiratCat P G).map φ) := by
    refine ⟨fun {W} u v huv => ?_⟩
    obtain ⟨Z, u₀, v₀, hu, hv⟩ := HomBirat.exists_rep_pair (P := P) (G := G) u v
    have hz : (toBiratCat P G).map Z.unop.hom.hom ≫ u = (toBiratCat P G).map u₀ := by
      rw [← hu]
      exact birat_toHom_comp_mk Z.unop.hom.hom Z.unop.hom.property.1
        Z.unop.hom.property.2 u₀
    have hz' : (toBiratCat P G).map Z.unop.hom.hom ≫ v = (toBiratCat P G).map v₀ := by
      rw [← hv]
      exact birat_toHom_comp_mk Z.unop.hom.hom Z.unop.hom.property.1
        Z.unop.hom.property.2 v₀
    have e2 : (toBiratCat P G).map (u₀ ≫ φ) = (toBiratCat P G).map (v₀ ≫ φ) := by
      rw [(toBiratCat P G).map_comp, (toBiratCat P G).map_comp, ← hz, ← hz']
      exact (Category.assoc _ _ _).trans
        ((congrArg (fun t => (toBiratCat P G).map Z.unop.hom.hom ≫ t) huv).trans
          (Category.assoc _ _ _).symm)
    have e3 : u₀ ≫ φ = v₀ ≫ φ := (toBiratCat P G).map_injective e2
    have e4 : u₀ = v₀ := (cancel_mono φ).mp e3
    rw [← hu, ← hv, e4]
  rw [hfeq]
  refine ⟨fun {W} u v huv => ?_⟩
  have h1 : (u ≫ a') ≫ (toBiratCat P G).map φ = (v ≫ a') ≫ (toBiratCat P G).map φ :=
    (Category.assoc _ _ _).trans (huv.trans (Category.assoc _ _ _).symm)
  have h2 : u ≫ a' = v ≫ a' := hmono2.right_cancellation _ _ h1
  have h3 : (u ≫ a') ≫ aa = (v ≫ a') ≫ aa := congrArg (fun t => t ≫ aa) h2
  have h4 : u ≫ (a' ≫ aa) = v ≫ (a' ≫ aa) :=
    (Category.assoc _ _ _).symm.trans (h3.trans (Category.assoc _ _ _))
  have h5 : u ≫ 𝟙 A = v ≫ 𝟙 A := by rw [← ha2]; exact h4
  exact ((Category.comp_id u).symm.trans h5).trans (Category.comp_id v)

/-! ## ★13. (iii)(b) —— co-angular pre-step があれば射は全部 co-angular

原文 (FrdI p.24):
> co-angular pre-step φ : A ->B, there exists a [uniquely determined] bijection of

★★★**鍵は「`𝒞^birat` の自己射はすべて co-angular」**である ——
標準形 `ψ = c' ≫ [ψ₀]` の添字の構造射 `c : A₃ ⟶ A` が
**`𝒞` の co-angular pre-step**で、`ψ₀` も `A₃ ⟶ A` だから、
★`𝒞` の (iii)(b) がそのまま当たる。 -/

include P in
/-- ★★★**`𝒞^birat` の自己射はすべて co-angular**。 -/
theorem birat_isCoAngular_endo {X : BiratCat P G} (ψ : X ⟶ X) :
    IsCoAngular (biratPre P G) ψ := by
  obtain ⟨A₃, c, ψ₀, cc, c', hcc, hcs, -, hc1, hc2, hψeq⟩ := birat_hom_repr P G ψ
  exact (birat_isCoAngular_repr P G ψ ψ₀ cc c' hc1 hc2 hψeq).mpr
    (G.core.coAngularOfPreStep c hcc hcs ψ₀)

include P in
/-- ★★★**[FrdI] Definition 1.3, (iii)(b)** の `𝒞^birat` 版。

★`α` は同型なので `φ = (φ ≫ α⁻¹) ≫ α` と書け、前半は自己射だから co-angular。 -/
theorem birat_coAngularOfPreStep {A' A : BiratCat P G} (α : A' ⟶ A)
    (hc : IsCoAngular (biratPre P G) α) (hs : IsPreStep (biratPre P G) α)
    (φ : A' ⟶ A) : IsCoAngular (biratPre P G) φ := by
  obtain ⟨p, hp1, hp2⟩ := birat_coaPre_inv P G α hc hs
  have heq : (φ ≫ p) ≫ α = φ :=
    (Category.assoc _ _ _).trans
      ((congrArg (fun t => φ ≫ t) hp2).trans (Category.comp_id φ))
  have hcomp : IsCoAngular (biratPre P G) ((φ ≫ p) ≫ α) :=
    birat_coAngularComp P G (φ ≫ p) α (birat_isCoAngular_endo P G (φ ≫ p)) hc
  rwa [heq] at hcomp

/-! ## ★14. isotropic 対象の辞書 —— (vii) を渡すための 3 補題

原文 (FrdI p.83):
> isotropic; ...

★★**`𝒞` の中で 2 つ、`𝒞^birat` との間で 1 つ**:

| 補題 | 内容 |
|---|---|
| `isCoAngular_of_isotropic` | isotropic 対象から出る pre-step は co-angular |
| `isotropic_of_coaPre` | ★**isotropic 対象への co-angular pre-step は、始域も isotropic** |
| `birat_isIsotropic_iff` | ★★`X` が `𝒞^birat` で isotropic ⟺ `𝒞` で isotropic |
-/

include G in
/-- ★★**isotropic 対象から出る pre-step は co-angular**。

★(v)(b) で `φ = β ≫ α`(`β` co-angular pre-step、`α` isometric pre-step)と分け、
`β` の終域も isotropic((vii)(b))なので `α` は同型になる。 -/
theorem isCoAngular_of_isotropic {A B : C} (hA : IsIsotropic P A) (φ : A ⟶ B)
    (hφ : IsPreStep P φ) : IsCoAngular P φ := by
  obtain ⟨M, β, α, hfac, hβc, hβs, hαi, hαs⟩ := G.core.preStepFactor φ hφ
  have hM : IsIsotropic P M := G.core.isotropicClosed β hA
  haveI hαiso : IsIso α := hM _ α hαi hαs
  rw [hfac]
  exact isCoAngular_compIso (Q := P) α (inv α) (IsIso.hom_inv_id α) (IsIso.inv_hom_id α) hβc

include G in
/-- ★★★**isotropic 対象への co-angular pre-step は、始域も isotropic**。

★★`A` の isotropic hull `ι : A ⟶ H` を取り、`B` が isotropic なので
`b = ι ≫ b̃` と分ける。★**`b` の co-angular 性を分解 `b = 𝟙 ≫ ι ≫ b̃` に当てる**と
`ι` が同型になる —— これが要点である。 -/
theorem isotropic_of_coaPre {A B : C} (b : A ⟶ B) (hbc : IsCoAngular P b)
    (hbs : IsPreStep P b) (hB : IsIsotropic P B) : IsIsotropic P A := by
  obtain ⟨H, ι, hιi, hιs, hιiso, hιuniv⟩ := G.core.isotropicHullExists A
  obtain ⟨bt, hbt, -⟩ := hιuniv B hB b
  have hdbt : P.degFr bt = 1 := by
    have h1 : P.degFr b = 1 := hbs.1
    rw [hbt, P.degFr_comp, show P.degFr ι = 1 from hιs.1, mul_one] at h1
    exact h1
  haveI hιIso : IsIso ι := by
    refine hbc A H (𝟙 A) ι bt ?_ hdbt hιi hιs
      (Or.inr (isBaseIsomorphism_of_isIso P (𝟙 A)))
    rw [Category.id_comp]
    exact hbt
  intro Dd h hisom hstep
  haveI : IsIso (inv ι ≫ h) :=
    hιiso Dd (inv ι ≫ h) ((isIsometric_of_isIso P (inv ι)).comp P hisom)
      ((isPreStep_of_isIso P (inv ι)).comp P hstep)
  have h4 : IsIso (ι ≫ (inv ι ≫ h)) := inferInstance
  rwa [← Category.assoc, IsIso.hom_inv_id, Category.id_comp] at h4

include P in
/-- ★★**`𝒞^birat` で isotropic なら `𝒞` でも isotropic**。

★`𝒞` の isometric pre-step `h` を `𝒞^birat` へ送ると同型になるので、
`h` は co-angular pre-step。★`Proposition 1.4, (iii)` で `𝒞` の中でも同型。 -/
theorem birat_isotropic_down {X : BiratCat P G}
    (hX : IsIsotropic (biratPre P G) X) : IsIsotropic P (biratDown P G X) := by
  intro W h hisom hstep
  haveI : IsIso ((toBiratCat P G).map h) :=
    hX (show BiratCat P G from W) ((toBiratCat P G).map h) (birat_isIsometric _)
      ((birat_isPreStep_iff (P := P) (G := G) h).mpr hstep)
  obtain ⟨hc, -⟩ := (birat_isIso_iff (P := P) (G := G) h).mp inferInstance
  exact prop_1_4_iii P G.core h ⟨hc, hisom⟩ hstep

include P in
/-- ★★★**`𝒞` で isotropic なら `𝒞^birat` でも isotropic**。

★標準形 `g = b' ≫ [g₀]` の `B₁` は `isotropic_of_coaPre` で isotropic、
そこから出る pre-step `g₀` は `isCoAngular_of_isotropic` で co-angular、
よって `[g₀]` は同型。 -/
theorem birat_isotropic_up {X : C} (hX : IsIsotropic P X) :
    IsIsotropic (biratPre P G) (show BiratCat P G from X) := by
  intro Y g _ hgs
  obtain ⟨B₁, b, g₀, bb, b', hbc, hbs, hbb, hb1, hb2, hgeq⟩ := birat_hom_repr P G g
  have hg₀s : IsPreStep P g₀ := by
    refine (birat_isPreStep_iff (P := P) (G := G) g₀).mp ?_
    have e1 : (toBiratCat P G).map g₀ = bb ≫ g := by
      rw [hgeq]
      exact (((Category.assoc _ _ _).symm.trans
        (congrArg (fun t => t ≫ (toBiratCat P G).map g₀) hb1)).trans
        (Category.id_comp _)).symm
    rw [e1]
    exact IsPreStep.comp (biratPre P G)
      (hbb ▸ (birat_isPreStep_iff (P := P) (G := G) b).mpr hbs) hgs
  have hB₁ : IsIsotropic P B₁ := isotropic_of_coaPre P G b hbc hbs hX
  have hg₀c : IsCoAngular P g₀ := isCoAngular_of_isotropic P G hB₁ g₀ hg₀s
  obtain ⟨mg, hmg⟩ : ∃ m : (show BiratCat P G from B₁) ⟶ Y,
      m = (toBiratCat P G).map g₀ := ⟨_, rfl⟩
  haveI : IsIso mg := hmg ▸ birat_isIso_of_coaPre g₀ hg₀c hg₀s
  haveI hb'iso : IsIso b' := ⟨⟨bb, hb2, hb1⟩⟩
  have hfin : IsIso (b' ≫ mg) := inferInstance
  rw [hgeq, ← hmg]
  exact hfin

include P in
/-- ★★★★**辞書 (iv) の「isotropic 対象」の条**。

原文 (FrdI p.83):
> isotropic; ...
-/
theorem birat_isIsotropic_iff (X : C) :
    IsIsotropic (biratPre P G) (show BiratCat P G from X) ↔ IsIsotropic P X :=
  ⟨fun h => birat_isotropic_down P G h, fun h => birat_isotropic_up P G h⟩

include P in
/-- ★★★**[FrdI] Definition 1.3, (vii)(b)** の `𝒞^birat` 版。

★辞書で `𝒞` へ落とし、標準形の `A₁` も isotropic にしてから
`𝒞` の (vii)(b) を当て、また辞書で戻す。 -/
theorem birat_isotropicClosed {A B : BiratCat P G} (φ : A ⟶ B)
    (hA : IsIsotropic (biratPre P G) A) : IsIsotropic (biratPre P G) B := by
  obtain ⟨A₁, a, φ₀, aa, a', hac, has, haa, ha1, ha2, hφeq⟩ := birat_hom_repr P G φ
  have hdA : IsIsotropic P (biratDown P G A) := birat_isotropic_down P G hA
  have hA₁ : IsIsotropic P A₁ := isotropic_of_coaPre P G a hac has hdA
  exact birat_isotropic_up P G (G.core.isotropicClosed φ₀ hA₁)

/-! ## ★15. (vii)(a) —— `𝒞^birat` の isotropic hull

原文 (FrdI p.25):
> (vii) (Isotropic Objects) (a) For every A ∈Ob(C), there exists an isotropic hull

★★**`𝒞` の isotropic hull `ι : A ⟶ H` をそのまま送る**。
`H` が `𝒞^birat` でも isotropic なのは `birat_isIsotropic_iff`。

★★★普遍性が要点である。標準形 `γ = a' ≫ [γ₀]` の `A₁` の hull `ι₁ : A₁ ⟶ H₁` を取り、
`a ≫ ι` を `ι₁ ≫ t` と分ける。★**`t` は isotropic 対象 `H₁` から出る pre-step なので
co-angular**(`isCoAngular_of_isotropic`)、よって `𝒞^birat` では**同型**であり、
`β := t⁻¹ ≫ [β₁]` が取れる。★一意性は `𝒞^birat` が totally epimorphic だから。 -/

include P in
/-- ★★★★**[FrdI] Definition 1.3, (vii)(a)** の `𝒞^birat` 版。 -/
theorem birat_isotropicHullExists (A : BiratCat P G) :
    ∃ (B : BiratCat P G) (φ : A ⟶ B), IsIsotropicHull (biratPre P G) φ := by
  obtain ⟨H, ι, hιi, hιs, hιiso, hιuniv⟩ := G.core.isotropicHullExists (biratDown P G A)
  refine ⟨(show BiratCat P G from H), (toBiratCat P G).map ι,
    birat_isIsometric _, (birat_isPreStep_iff (P := P) (G := G) ι).mpr hιs,
    birat_isotropic_up P G hιiso, ?_⟩
  intro Cc hCc γ
  have hdCc : IsIsotropic P (biratDown P G Cc) := birat_isotropic_down P G hCc
  obtain ⟨A₁, a, γ₀, aa, a', hac, has, haa, ha1, ha2, hγeq⟩ := birat_hom_repr P G γ
  obtain ⟨H₁, ι₁, hι₁i, hι₁s, hι₁iso, hι₁univ⟩ := G.core.isotropicHullExists A₁
  obtain ⟨β₁, hβ₁, -⟩ := hι₁univ (biratDown P G Cc) hdCc γ₀
  obtain ⟨t, ht, -⟩ := hι₁univ H hιiso (a ≫ ι)
  -- ★`t` は pre-step
  have hdt : P.degFr t = 1 := by
    have h1 : P.degFr (a ≫ ι) = P.degFr (ι₁ ≫ t) := congrArg P.degFr ht
    rw [P.degFr_comp, P.degFr_comp, show P.degFr a = 1 from has.1,
      show P.degFr ι = 1 from hιs.1, show P.degFr ι₁ = 1 from hι₁s.1] at h1
    simpa using h1.symm
  have hbt : IsBaseIsomorphism P t := by
    have hb1 : IsIso (P.Base (ι₁ ≫ t)) :=
      ht ▸ isBaseIsomorphism_comp P has.2 hιs.2
    haveI : IsIso (P.Base ι₁) := hι₁s.2
    haveI : IsIso (P.Base ι₁ ≫ P.Base t) := by rwa [P.Base_comp] at hb1
    exact IsIso.of_isIso_comp_left (P.Base ι₁) (P.Base t)
  have hts : IsPreStep P t := ⟨hdt, hbt⟩
  -- ★★`t` は isotropic 対象から出る pre-step なので co-angular、よって `𝒞^birat` で同型
  have htc : IsCoAngular P t := isCoAngular_of_isotropic P G hι₁iso t hts
  obtain ⟨mt, hmt⟩ : ∃ m : (show BiratCat P G from H₁) ⟶ (show BiratCat P G from H),
      m = (toBiratCat P G).map t := ⟨_, rfl⟩
  haveI hmtiso : IsIso mt := hmt ▸ birat_isIso_of_coaPre t htc hts
  obtain ⟨mt', hmt1, hmt2⟩ :
      ∃ m' : (show BiratCat P G from H) ⟶ (show BiratCat P G from H₁),
        mt ≫ m' = 𝟙 _ ∧ m' ≫ mt = 𝟙 _ := hmtiso.out
  obtain ⟨mβ, hmβ⟩ : ∃ m : (show BiratCat P G from H₁) ⟶ Cc,
      m = (toBiratCat P G).map β₁ := ⟨_, rfl⟩
  obtain ⟨mι, hmι⟩ : ∃ m : A ⟶ (show BiratCat P G from H),
      m = (toBiratCat P G).map ι := ⟨_, rfl⟩
  obtain ⟨mι₁, hmι₁⟩ : ∃ m : (show BiratCat P G from A₁) ⟶ (show BiratCat P G from H₁),
      m = (toBiratCat P G).map ι₁ := ⟨_, rfl⟩
  obtain ⟨mγ₀, hmγ₀⟩ : ∃ m : (show BiratCat P G from A₁) ⟶ Cc,
      m = (toBiratCat P G).map γ₀ := ⟨_, rfl⟩
  have mc1 : (toBiratCat P G).map (ι₁ ≫ t)
      = (toBiratCat P G).map ι₁ ≫ (toBiratCat P G).map t := (toBiratCat P G).map_comp _ _
  have mc2 : (toBiratCat P G).map (a ≫ ι)
      = (toBiratCat P G).map a ≫ (toBiratCat P G).map ι := (toBiratCat P G).map_comp _ _
  have mc3 : (toBiratCat P G).map (ι₁ ≫ β₁)
      = (toBiratCat P G).map ι₁ ≫ (toBiratCat P G).map β₁ := (toBiratCat P G).map_comp _ _
  have e_t : mι₁ ≫ mt = aa ≫ mι := by
    rw [hmι₁, hmt, hmι, haa]
    exact (mc1.symm.trans (congrArg (toBiratCat P G).map ht.symm)).trans mc2
  have e_γ : mι₁ ≫ mβ = mγ₀ := by
    rw [hmι₁, hmβ, hmγ₀]
    exact mc3.symm.trans (congrArg (toBiratCat P G).map hβ₁.symm)
  have hγ2 : γ = a' ≫ mγ₀ := by rw [hmγ₀]; exact hγeq
  -- ★等式の組み立て
  have q2 : (aa ≫ mι) ≫ mt' = mι₁ := by
    have r1 : (mι₁ ≫ mt) ≫ mt' = mι₁ ≫ (mt ≫ mt') := Category.assoc _ _ _
    have r2 : mι₁ ≫ (mt ≫ mt') = mι₁ :=
      (congrArg (fun z => mι₁ ≫ z) hmt1).trans (Category.comp_id _)
    exact ((congrArg (fun z => z ≫ mt') e_t).symm.trans r1).trans r2
  have s0 : (aa ≫ mι) ≫ (mt' ≫ mβ) = mι₁ ≫ mβ :=
    (Category.assoc _ _ _).symm.trans (congrArg (fun z => z ≫ mβ) q2)
  have s1 : aa ≫ (mι ≫ (mt' ≫ mβ)) = mγ₀ :=
    ((Category.assoc _ _ _).symm.trans s0).trans e_γ
  have s2 : γ = mι ≫ (mt' ≫ mβ) := by
    have u1 : a' ≫ mγ₀ = a' ≫ (aa ≫ (mι ≫ (mt' ≫ mβ))) :=
      congrArg (fun z => a' ≫ z) s1.symm
    have u2 : a' ≫ (aa ≫ (mι ≫ (mt' ≫ mβ))) = (a' ≫ aa) ≫ (mι ≫ (mt' ≫ mβ)) :=
      (Category.assoc _ _ _).symm
    have u3 : (a' ≫ aa) ≫ (mι ≫ (mt' ≫ mβ)) = 𝟙 A ≫ (mι ≫ (mt' ≫ mβ)) :=
      congrArg (fun z => z ≫ (mι ≫ (mt' ≫ mβ))) ha2
    exact ((hγ2.trans u1).trans (u2.trans u3)).trans (Category.id_comp _)
  have hs2' : γ = (toBiratCat P G).map ι ≫ (mt' ≫ mβ) := by rw [← hmι]; exact s2
  refine ⟨mt' ≫ mβ, hs2', ?_⟩
  intro y hy
  have hy' : γ = (toBiratCat P G).map ι ≫ y := hy
  haveI hepi : Epi ((toBiratCat P G).map ι) :=
    birat_totEpi P G _ _ ((toBiratCat P G).map ι)
  exact hepi.left_cancellation _ _ (hy'.symm.trans hs2')

end ABC3.Found.FrdI
