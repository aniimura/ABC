import ABC3.Found.FrdI.Def31
import Mathlib.CategoryTheory.Conj

/-!
# [FrdI] Theorem 3.4 —— 底と Frobenius 次数の圏論性

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.62–p.69。

原文 (FrdI p.62):
> (Category-theoreticity of the Base and Frobenius De-

原文 (FrdI p.63):
> First, we consider assertion (i). Since iso-subanchors are manifestly pre-

## ★この定理の規模(測定)

**5 条**あり、主張は多い。★この定理は **IUTchI・IUTchII が直接引く**もので、
`Remark 3.4.1` と `Proposition 3.11, (ii)(iii)` の入口でもある。

## ★★(i) の第 1 歩 —— ここで実装する部分

原文 (FrdI p.63):
> served by any equivalence of categories, it follows from our assumption that C

★原文は「**iso-subanchor は任意の圏同値で明らかに保たれる**」と一言で済ませる。
★★**「明らかに」の中身は 5 段**である:

| 段 | 内容 |
|---|---|
| 1 | **irreducible 射**が保たれる |
| 2 | **anchor** が保たれる |
| 3 | **subanchor** が保たれる |
| 4 | **categorical quotient** と **mono-minimal** が保たれる |
| 5 | ゆえに **iso-subanchor** が保たれ、`quasi-isotropic 型`から isotropic 対象が保たれる |

★段 1–3 をここで取る。★★**要点は「圏同値は分解を反射する」**こと ——
充満忠実性で射を戻し、本質的全射性で中間対象を戻す。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v₁ u₁ v₂ u₂

variable {C : Type u₁} [Category.{v₁} C] {D : Type u₂} [Category.{v₂} D]
  (Ψ : C ⥤ D) [Ψ.IsEquivalence]

/-! ## ★段 1 —— irreducible 射は保たれる -/

/-- ★★**圏同値は分解を反射する** —— `Ψ.map φ = β ≫ α` なら、
`φ = β₀ ≫ α₀` と分解でき、`β`・`α` はそれぞれ `β₀`・`α₀` の像に(同型を除いて)一致する。

★**本質的全射性**で中間対象 `X` を `Ψ.obj X₀` に戻し、
★**充満性**で `β`・`α` を戻し、★**忠実性**で分解の等式を `C` の側へ落とす。 -/
theorem exists_factor_of_map_factor {A B : C} (φ : A ⟶ B) {X : D}
    (β : Ψ.obj A ⟶ X) (α : X ⟶ Ψ.obj B) (h : Ψ.map φ = β ≫ α) :
    ∃ (X₀ : C) (β₀ : A ⟶ X₀) (α₀ : X₀ ⟶ B) (e : Ψ.obj X₀ ≅ X),
      φ = β₀ ≫ α₀ ∧ β = Ψ.map β₀ ≫ e.hom ∧ α = e.inv ≫ Ψ.map α₀ := by
  obtain ⟨X₀, ⟨e⟩⟩ : ∃ X₀ : C, Nonempty (Ψ.obj X₀ ≅ X) := ⟨Ψ.objPreimage X, ⟨Ψ.objObjPreimageIso X⟩⟩
  obtain ⟨β₀, hβ₀⟩ := Ψ.map_surjective (β ≫ e.inv)
  obtain ⟨α₀, hα₀⟩ := Ψ.map_surjective (e.hom ≫ α)
  refine ⟨X₀, β₀, α₀, e, ?_, ?_, ?_⟩
  · refine Ψ.map_injective ?_
    rw [Ψ.map_comp, hβ₀, hα₀, h]
    simp
  · rw [hβ₀]; simp
  · rw [hα₀]; simp

/-- ★★**圏同値は irreducible 射を保つ**。

★`¬ IsIso` は充満忠実性(同型を反射する)から、
★分解の条件は `exists_factor_of_map_factor` から。 -/
theorem isIrreducibleMor_map {A B : C} {φ : A ⟶ B} (h : IsIrreducibleMor φ) :
    IsIrreducibleMor (Ψ.map φ) := by
  refine ⟨fun hiso => h.1 ?_, ?_⟩
  · haveI := hiso
    exact (Functor.FullyFaithful.ofFullyFaithful Ψ).isIso_of_isIso_map φ
  · intro X β α hfac
    obtain ⟨X₀, β₀, α₀, e, hfac₀, hβ, hα⟩ :=
      exists_factor_of_map_factor Ψ φ β α hfac
    rcases h.2 X₀ β₀ α₀ hfac₀ with hb | ha
    · left
      haveI := hb
      rw [hβ]
      infer_instance
    · right
      haveI := ha
      rw [hα]
      infer_instance

/-! ## ★段 2・3 —— anchor と subanchor は保たれる -/

/-- ★★**圏同値は anchor を保つ**。

★`A` から出る irreducible 射の同型類が有限個であること —— これは
`Under A` の有限集合で押さえられている。★**`Ψ` で写した先でも、
戻して押さえ、また写せばよい**。 -/
theorem isAnchor_map {A : C} (h : IsAnchor C A) : IsAnchor D (Ψ.obj A) := by
  obtain ⟨s, hsfin, hs⟩ := h
  refine ⟨(fun Z : Under A => Under.mk (Ψ.map Z.hom)) '' s, hsfin.image _, ?_⟩
  intro B φ hφ
  -- ★`φ` を `C` へ戻す
  obtain ⟨B₀, ⟨e⟩⟩ : ∃ B₀ : C, Nonempty (Ψ.obj B₀ ≅ B) :=
    ⟨Ψ.objPreimage B, ⟨Ψ.objObjPreimageIso B⟩⟩
  obtain ⟨φ₀, hφ₀⟩ := Ψ.map_surjective (φ ≫ e.inv)
  have hφ₀irr : IsIrreducibleMor φ₀ := by
    refine ⟨fun hiso => hφ.1 ?_, ?_⟩
    · haveI : IsIso (Ψ.map φ₀) := by haveI := hiso; infer_instance
      have : IsIso (φ ≫ e.inv) := hφ₀ ▸ inferInstance
      haveI := this
      have h2 : φ = (φ ≫ e.inv) ≫ e.hom := by simp
      rw [h2]
      infer_instance
    · intro X₀ β₀ α₀ hfac
      refine hφ.2 (Ψ.obj X₀) (Ψ.map β₀) (Ψ.map α₀ ≫ e.hom) ?_ |>.imp ?_ ?_
      · rw [← Category.assoc, ← Ψ.map_comp, ← hfac, hφ₀]
        simp
      · intro hb
        exact (Functor.FullyFaithful.ofFullyFaithful Ψ).isIso_of_isIso_map β₀
      · intro ha
        haveI := ha
        haveI : IsIso (Ψ.map α₀ ≫ e.hom) := ha
        haveI : IsIso (Ψ.map α₀) := IsIso.of_isIso_comp_right _ e.hom
        exact (Functor.FullyFaithful.ofFullyFaithful Ψ).isIso_of_isIso_map α₀
  obtain ⟨Z, hZs, ⟨j⟩⟩ := hs B₀ φ₀ hφ₀irr
  refine ⟨Under.mk (Ψ.map Z.hom), ⟨Z, hZs, rfl⟩, ⟨?_⟩⟩
  -- `Under.mk φ ≅ Under.mk (Ψ.map Z.hom)`
  refine Under.isoMk (e.symm ≪≫ (Ψ.mapIso ((Under.forget A).mapIso j))) ?_
  have hw : φ₀ ≫ (Under.forget A).map j.hom = Z.hom :=
    Under.w j.hom
  show φ ≫ (e.inv ≫ Ψ.map ((Under.forget A).map j.hom)) = Ψ.map Z.hom
  rw [← Category.assoc, ← hφ₀, ← Ψ.map_comp, hw]
  rfl


/-! ## ★段 4 —— categorical quotient と mono-minimal は保たれる -/

/-- ★★**圏同値は `Aut` を全単射に写す**。

★単射性は忠実性、★全射性は充満性(`preimageIso`)から。
★★**mono-minimal の移送で部分群を戻すときに要る**。 -/
theorem mapAut_bijective (A : C) : Function.Bijective (Ψ.mapAut A) := by
  constructor
  · intro g₁ g₂ hg
    refine Aut.ext ?_
    refine Ψ.map_injective ?_
    exact congrArg (fun z : Aut (Ψ.obj A) => z.hom) hg
  · intro g
    refine ⟨(Functor.FullyFaithful.ofFullyFaithful Ψ).preimageIso g, ?_⟩
    refine Aut.ext ?_
    show Ψ.map ((Functor.FullyFaithful.ofFullyFaithful Ψ).preimageIso g).hom = g.hom
    simp

/-- ★**同型として述べたもの** —— `Aut A ≃* Aut (Ψ.obj A)`。 -/
noncomputable def mapAutEquiv (A : C) : Aut A ≃* Aut (Ψ.obj A) :=
  MulEquiv.ofBijective (Ψ.mapAut A) (mapAut_bijective Ψ A)

/-- ★★**圏同値は categorical quotient を保つ**。

★不変性は `Ψ.map` の関手性から。★**普遍性の側**が要点で、
`Cc` を `Ψ.obj Cc₀` に戻し(本質的全射性)、`ψ` を戻し(充満性)、
`C` の側の一意性を使って戻す。 -/
theorem isCategoricalQuotient_map {A B : C} (G : Subgroup (Aut A)) {φ : A ⟶ B}
    (h : IsCategoricalQuotient G φ) :
    IsCategoricalQuotient (G.map (Ψ.mapAut A)) (Ψ.map φ) := by
  constructor
  · rintro γ' ⟨γ, hγ, rfl⟩
    show (Ψ.mapAut A γ : Aut (Ψ.obj A)).hom ≫ Ψ.map φ = Ψ.map φ
    show Ψ.map (γ : Aut A).hom ≫ Ψ.map φ = Ψ.map φ
    rw [← Ψ.map_comp, h.1 γ hγ]
  · intro Cc ψ hψ
    obtain ⟨Cc₀, ⟨e⟩⟩ : ∃ Cc₀ : C, Nonempty (Ψ.obj Cc₀ ≅ Cc) :=
      ⟨Ψ.objPreimage Cc, ⟨Ψ.objObjPreimageIso Cc⟩⟩
    obtain ⟨ψ₀, hψ₀⟩ := Ψ.map_surjective (ψ ≫ e.inv)
    have hψ₀inv : ∀ γ ∈ G, ((γ : Aut A).hom : A ⟶ A) ≫ ψ₀ = ψ₀ := by
      intro γ hγ
      refine Ψ.map_injective ?_
      rw [Ψ.map_comp, hψ₀]
      have hthis : Ψ.map ((γ : Aut A).hom : A ⟶ A) ≫ ψ = ψ :=
        hψ (Ψ.mapAut A γ) ⟨γ, hγ, rfl⟩
      show Ψ.map (γ : Aut A).hom ≫ (ψ ≫ e.inv) = ψ ≫ e.inv
      rw [← Category.assoc, hthis]
    obtain ⟨ψ₀', hψ₀'eq, hψ₀'uniq⟩ := h.2 Cc₀ ψ₀ hψ₀inv
    refine ⟨Ψ.map ψ₀' ≫ e.hom, ?_, ?_⟩
    · show ψ = Ψ.map φ ≫ (Ψ.map ψ₀' ≫ e.hom)
      rw [← Category.assoc, ← Ψ.map_comp, ← hψ₀'eq, hψ₀]
      simp
    · intro y hy
      have hy₀ : ψ₀ = φ ≫ (Ψ.preimage (y ≫ e.inv)) := by
        refine Ψ.map_injective ?_
        rw [Ψ.map_comp, Ψ.map_preimage, hψ₀, hy]
        simp
      have := hψ₀'uniq _ hy₀
      rw [← this]
      simp

/-- ★**`Aut` の移送** —— `Ψ` による写しと、同型による共役の合成。 -/
noncomputable def autTransfer {A₀ : C} {A'' : D} (ε : Ψ.obj A₀ ≅ A'') :
    Aut A₀ ≃* Aut A'' :=
  (mapAutEquiv Ψ A₀).trans ε.conjAut

/-- ★★★**圏同値は mono-minimal categorical quotient を保つ**。

★★**要点は `G'` と `e` が定義で全称量化されている**こと ——
移送では **`C` 側の部分群 `G₀` と同型 `e₀` をこちらが構成してよい**。 -/
theorem isMonoMinimalQuotient_map {A B : C} (G : Subgroup (Aut A)) {φ : A ⟶ B}
    (h : IsMonoMinimalQuotient G φ) :
    IsMonoMinimalQuotient (G.map (Ψ.mapAut A)) (Ψ.map φ) := by
  intro A'' ζ φ'' hmono hfac G'' e hcompat
  -- ★手 1: 対象と射を `C` へ戻す
  set A₀ := Ψ.objPreimage A'' with hA₀
  set ε : Ψ.obj A₀ ≅ A'' := Ψ.objObjPreimageIso A'' with hε
  set ζ₀ : A ⟶ A₀ := Ψ.preimage (ζ ≫ ε.inv) with hζ₀def
  have hζ₀ : Ψ.map ζ₀ = ζ ≫ ε.inv := Ψ.map_preimage _
  haveI : Mono (ζ ≫ ε.inv) := by haveI := hmono; infer_instance
  have hmono₀ : Mono ζ₀ := by
    haveI : Mono (Ψ.map ζ₀) := hζ₀ ▸ inferInstance
    exact Ψ.mono_of_mono_map this
  set φ₀ : A₀ ⟶ B := Ψ.preimage (ε.hom ≫ φ'') with hφ₀def
  have hfac₀ : φ = ζ₀ ≫ φ₀ := by
    refine Ψ.map_injective ?_
    rw [Ψ.map_comp, hζ₀, hφ₀def, Ψ.map_preimage, hfac]
    simp
  -- ★手 2: 群を戻す
  set K : Aut A₀ ≃* Aut A'' := autTransfer Ψ ε with hK
  set i₁ : G ≃* G.map (Ψ.mapAut A) :=
    Subgroup.equivMapOfInjective G (Ψ.mapAut A) (mapAut_bijective Ψ A).1 with hi₁
  set G₀ : Subgroup (Aut A₀) := G''.map K.symm.toMonoidHom with hG₀
  set e₀ : G ≃* G₀ := i₁.trans (e.trans (K.symm.subgroupMap G'')) with he₀
  -- ★手 3: 両立条件
  have hcompat₀ : ∀ γ : G, ((γ : Aut A).hom : A ⟶ A) ≫ ζ₀
      = ζ₀ ≫ ((e₀ γ : Aut A₀).hom : A₀ ⟶ A₀) := by
    intro γ
    refine Ψ.map_injective ?_
    have hKe : (K (e₀ γ : Aut A₀) : Aut A'').hom = ((e (i₁ γ) : Aut A'')).hom :=
      congrArg (fun z : Aut A'' => z.hom) (K.apply_symm_apply _)
    have hKf : (K (e₀ γ : Aut A₀) : Aut A'').hom
        = ε.inv ≫ Ψ.map ((e₀ γ : Aut A₀).hom) ≫ ε.hom := rfl
    have hval : Ψ.map ((e₀ γ : Aut A₀).hom)
        = ε.hom ≫ ((e (i₁ γ) : Aut A'').hom) ≫ ε.inv := by
      rw [← hKe, hKf]
      simp
    have hc : Ψ.map ((γ : Aut A).hom : A ⟶ A) ≫ ζ
        = ζ ≫ ((e (i₁ γ) : Aut A'').hom) := hcompat (i₁ γ)
    rw [Ψ.map_comp, Ψ.map_comp, hζ₀, hval, ← Category.assoc, hc]
    simp
  -- ★手 4: `C` の側で結論を得て戻す
  haveI : IsIso ζ₀ := h A₀ ζ₀ φ₀ hmono₀ hfac₀ G₀ e₀ hcompat₀
  haveI : IsIso (Ψ.map ζ₀) := inferInstance
  haveI : IsIso (ζ ≫ ε.inv) := hζ₀ ▸ inferInstance
  have hz : ζ = (ζ ≫ ε.inv) ≫ ε.hom := by simp
  rw [hz]
  infer_instance

/-! ### ★★★残り —— isotropic 対象が保たれること

★★上で **iso-subanchor が保たれる**ことが取れた。★`quasi-isotropic 型`の定義は

  `∀ A, ¬ IsIsotropic P A ↔ IsIsoSubanchor C A`

なので、両側でこれを使えば **`Ψ` は isotropic 対象を保つ**が出る。

★**残る部品は 1 つ**: 逆向き(`IsIsoSubanchor D (Ψ.obj A) → IsIsoSubanchor C A`)。
準逆関手 `Ψ.inv` に上の補題を当てると
`IsIsoSubanchor C (Ψ.inv.obj (Ψ.obj A))` が出るので、
★**`IsIsoSubanchor` が同型で不変**であること(`φ' ≫ iso.inv` を取る)を
足せばよい。categorical quotient と mono-minimal が
**同型との後合成で保たれる**ことを見ればよい。
-/

/-! ### ★★★段 5 の設計(実装済み。以下は記録)

★残るのは **mono-minimal の移送**である。定義

  `IsMonoMinimalQuotient G φ := ∀ A' ζ φ', Mono ζ → φ = ζ ≫ φ' →`
  `  ∀ G' (e : G ≃* G'), (両立) → IsIso ζ`

で **`G'` と `e` が全称量化**されているので、移送では **`C` 側の `G₀`・`e₀` を
こちらが構成してよい**。これが要点である。

**手順**:
1. `A'' ≅ Ψ.obj A₀`(本質的全射性)、`ζ₀ := Ψ.preimage (ζ ≫ ε.inv)`、
   `φ₀ := Ψ.preimage (ε.hom ≫ φ'')`。★`Mono ζ₀` は**忠実関手が mono を反射する**
   ことから(`Ψ.map ζ₀ = ζ ≫ ε.inv` は mono ≫ 同型)
2. `K : Aut A₀ ≃* Aut A''` を **`Ψ.mapAut A₀` と `ε.conjAut` の合成**として作る
   (どちらも全単射 —— 前者は `mapAut_bijective`、後者は `Iso.conjAut`)
3. `G₀ := G''.comap K`、`e₀ : G ≃* G₀` は
   `G ≃* G.map (Ψ.mapAut A)`(単射)→ `e` → `K.symm` の制限、で作る
4. 両立条件は `Ψ.map` で写して確かめる。★計算すると
   `Ψ.map (e₀ γ).hom = ε.hom ≫ (e γ').hom ≫ ε.inv` なので、
   両辺とも `ζ ≫ (e γ').hom ≫ ε.inv` に落ちる
5. `h` を当てて `IsIso ζ₀` ⟹ `IsIso (Ψ.map ζ₀) = IsIso (ζ ≫ ε.inv)` ⟹ `IsIso ζ`

★★これが済めば `IsIsoSubanchor` の移送が
`isSubanchor_map` ＋ `isCategoricalQuotient_map` ＋ 本補題で組み上がり、
★★★**`quasi-isotropic 型`の定義(`¬IsIsotropic ↔ IsIsoSubanchor`)から
「`Ψ` は isotropic 対象を保つ」**——`Theorem 3.4, (i)` の第 1 歩——が出る。
-/

/-- ★**圏同値は subanchor を保つ**。 -/
theorem isSubanchor_map {A : C} (h : IsSubanchor C A) : IsSubanchor D (Ψ.obj A) := by
  obtain ⟨B, φ, hB⟩ := h
  exact ⟨Ψ.obj B, Ψ.map φ, isAnchor_map Ψ hB⟩

/-- ★★★★★**圏同値は iso-subanchor を保つ** —— 原文が
「iso-subanchors are manifestly preserved by any equivalence of categories」と
一言で済ませる主張。

★★**「明らかに」の中身は 5 段**だった:
irreducible 射 → anchor → subanchor → categorical quotient → mono-minimal。 -/
theorem isIsoSubanchor_map {A : C} (h : IsIsoSubanchor C A) :
    IsIsoSubanchor D (Ψ.obj A) := by
  obtain ⟨B, hB, G, φ, hcq, hmm⟩ := h
  exact ⟨Ψ.obj B, isSubanchor_map Ψ hB, G.map (Ψ.mapAut B), Ψ.map φ,
    isCategoricalQuotient_map Ψ G hcq, isMonoMinimalQuotient_map Ψ G hmm⟩


/-! ## ★段 5 の仕上げ —— isotropic 対象が保たれる -/

omit [Ψ.IsEquivalence] in
/-- ★**`IsIsoSubanchor` は同型で不変** —— `φ ≫ k.inv` を取ればよい。 -/
theorem isIsoSubanchor_of_iso {A A' : C} (k : A ≅ A') (h : IsIsoSubanchor C A') :
    IsIsoSubanchor C A := by
  obtain ⟨B, hB, G, φ, hcq, hmm⟩ := h
  refine ⟨B, hB, G, φ ≫ k.inv, ⟨?_, ?_⟩, ?_⟩
  · intro γ hγ
    rw [← Category.assoc, hcq.1 γ hγ]
  · intro Cc ψ hψ
    obtain ⟨ψ'', hψ''eq, hψ''uniq⟩ := hcq.2 Cc ψ hψ
    refine ⟨k.hom ≫ ψ'', ?_, ?_⟩
    · show ψ = (φ ≫ k.inv) ≫ (k.hom ≫ ψ'')
      rw [Category.assoc, ← Category.assoc k.inv, k.inv_hom_id, Category.id_comp]
      exact hψ''eq
    · intro y hy
      have hy2 : ψ = φ ≫ (k.inv ≫ y) := by
        rw [← Category.assoc]; exact hy
      have hu := hψ''uniq (k.inv ≫ y) hy2
      rw [← hu]
      simp
  · intro A'' ζ φ'' hmono hfac G' e hcompat
    refine hmm A'' ζ (φ'' ≫ k.hom) hmono ?_ G' e hcompat
    rw [← Category.assoc, ← hfac, Category.assoc, k.inv_hom_id, Category.comp_id]

/-- ★★**`IsIsoSubanchor` は圏同値で「保たれ、かつ反射される」**。 -/
theorem isIsoSubanchor_map_iff {A : C} :
    IsIsoSubanchor C A ↔ IsIsoSubanchor D (Ψ.obj A) := by
  refine ⟨isIsoSubanchor_map Ψ, fun h => ?_⟩
  have h' : IsIsoSubanchor C (Ψ.inv.obj (Ψ.obj A)) := isIsoSubanchor_map Ψ.inv h
  exact isIsoSubanchor_of_iso ((Ψ.asEquivalence.unitIso).app A) h'

section Isotropic

universe v₃ u₃ w₃ v₄ u₄ w₄

variable {E₁ : Type u₃} [Category.{v₃} E₁] {Φ₁ : MonoidOn.{v₃, u₃, w₃} E₁}
  (P₁ : PreFrobenioid C Φ₁)
  {E₂ : Type u₄} [Category.{v₄} E₂] {Φ₂ : MonoidOn.{v₄, u₄, w₄} E₂}
  (P₂ : PreFrobenioid D Φ₂)

/-- ★★★★★**[FrdI] Theorem 3.4, (i) の第 1 主張** —— `𝒞₁`・`𝒞₂` が
quasi-isotropic 型なら、**圏同値 `Ψ` は isotropic 対象を保つ**(かつ反射する)。

原文 (FrdI p.63):
> served by any equivalence of categories, it follows from our assumption that C

★★`quasi-isotropic 型`の定義は `¬ IsIsotropic ↔ IsIsoSubanchor` であり、
★**`IsIsoSubanchor` は純粋に圏論的**なので圏同値で保たれる
(`isIsoSubanchor_map_iff`)。それだけで出る。 -/
theorem thm_3_4_i_isotropic (h₁ : IsOfQuasiIsotropicType C P₁)
    (h₂ : IsOfQuasiIsotropicType D P₂) (A : C) :
    IsIsotropic P₁ A ↔ IsIsotropic P₂ (Ψ.obj A) := by
  rw [← not_iff_not, h₁ A, h₂ (Ψ.obj A)]
  exact isIsoSubanchor_map_iff Ψ

/-- ★★★★**[FrdI] Theorem 3.4, (i) の第 2 主張** —— `Ψ` は **isotropic hull を保つ**。

原文 (FrdI p.63):
> mainder of assertion (i) follows formally from [the definitions and] Proposition

★★鍵は `Proposition 1.9, (vi)` の**圏論的な特徴づけ**

  `IsIsotropicHull φ ↔ (終域が isotropic) ∧ (isotropic 始域の射に関して minimal-coadjoint)`

である。★**右辺は「isotropic 対象」という語彙だけで書けている**ので、
第 1 主張(`thm_3_4_i_isotropic`)がそのまま効く。 -/
theorem isotropicHull_map (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (h₁ : IsOfQuasiIsotropicType C P₁) (h₂ : IsOfQuasiIsotropicType D P₂)
    {A B : C} {φ : A ⟶ B} (h : IsIsotropicHull P₁ φ) :
    IsIsotropicHull P₂ (Ψ.map φ) := by
  rw [prop_1_9_vi P₂ F₂]
  rw [prop_1_9_vi P₁ F₁] at h
  refine ⟨(thm_3_4_i_isotropic Ψ P₁ P₂ h₁ h₂ B).mp h.1, ?_⟩
  intro X α β hfac hX
  obtain ⟨X₀, α₀, β₀, e, hfac₀, hα, hβ⟩ := exists_factor_of_map_factor Ψ φ α β hfac
  have hX₀ : IsIsotropic P₂ (Ψ.obj X₀) := isIsotropic_of_iso P₂ e hX
  have hX₀' : IsIsotropic P₁ X₀ := (thm_3_4_i_isotropic Ψ P₁ P₂ h₁ h₂ X₀).mpr hX₀
  haveI : IsIso β₀ := h.2 X₀ α₀ β₀ hfac₀ hX₀'
  rw [hβ]
  infer_instance

/-- ★★★★**[FrdI] Theorem 3.4, (i) の第 3 主張** —— `Ψ` は **等長 pre-step を保つ**。

★★鍵は `Proposition 1.9, (vii)` の**圏論的な特徴づけ**

  `φ` が等長 pre-step ⟺ `φ ≫ (B の hull)` が `A` の hull

である。★hull の保存(第 2 主張)から直ちに従う。 -/
theorem isometricPreStep_map (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (h₁ : IsOfQuasiIsotropicType C P₁) (h₂ : IsOfQuasiIsotropicType D P₂)
    {A B : C} {φ : A ⟶ B} (hφi : IsIsometric P₁ φ) (hφs : IsPreStep P₁ φ) :
    IsIsometric P₂ (Ψ.map φ) ∧ IsPreStep P₂ (Ψ.map φ) := by
  obtain ⟨Bi, ψ, hψ⟩ := F₁.isotropicHullExists B
  have hcomp : IsIsotropicHull P₁ (φ ≫ ψ) :=
    (prop_1_9_vii P₁ F₁ φ ψ hψ).mp ⟨hφi, hφs⟩
  have hmapψ : IsIsotropicHull P₂ (Ψ.map ψ) :=
    isotropicHull_map Ψ P₁ P₂ F₁ F₂ h₁ h₂ hψ
  have hmapcomp : IsIsotropicHull P₂ (Ψ.map φ ≫ Ψ.map ψ) := by
    rw [← Ψ.map_comp]
    exact isotropicHull_map Ψ P₁ P₂ F₁ F₂ h₁ h₂ hcomp
  exact (prop_1_9_vii P₂ F₂ (Ψ.map φ) (Ψ.map ψ) hmapψ).mpr hmapcomp

/-- ★★★★**[FrdI] Theorem 3.4, (i) の第 4 主張(関手の存在)** ——
`Ψ` は **`Ψ^istr : 𝒞₁^istr ⥤ 𝒞₂^istr` に制限される**。

原文 (FrdI p.62):
> there exists a 1-unique functor Ψistr : Cistr

★★対象の側は第 1 主張(`Ψ` は isotropic 対象を保つ)、
射の側は `Ψ.map` そのもの。★**充満部分圏なので関手性は自明に落ちる。** -/
def psiIstr (h₁ : IsOfQuasiIsotropicType C P₁) (h₂ : IsOfQuasiIsotropicType D P₂) :
    Istr P₁ ⥤ Istr P₂ where
  obj A := ⟨Ψ.obj A.obj, (thm_3_4_i_isotropic Ψ P₁ P₂ h₁ h₂ A.obj).mp A.property⟩
  map {_ _} f := ObjectProperty.homMk (Ψ.map f.hom)
  map_id _ := InducedCategory.hom_ext (Ψ.map_id _)
  map_comp _ _ := InducedCategory.hom_ext (Ψ.map_comp _ _)

/-- ★**`Ψ^istr` は忠実**(`Ψ` が忠実だから)。 -/
theorem psiIstr_faithful (h₁ : IsOfQuasiIsotropicType C P₁)
    (h₂ : IsOfQuasiIsotropicType D P₂) :
    (psiIstr Ψ P₁ P₂ h₁ h₂).Faithful where
  map_injective {_ _ f g} h := by
    refine InducedCategory.hom_ext (Ψ.map_injective ?_)
    exact congrArg (fun m : (psiIstr Ψ P₁ P₂ h₁ h₂).obj _ ⟶ _ => m.hom) h

/-- ★**`Ψ^istr` は充満**(`Ψ` が充満で、`𝒞^istr` が充満部分圏だから)。 -/
theorem psiIstr_full (h₁ : IsOfQuasiIsotropicType C P₁)
    (h₂ : IsOfQuasiIsotropicType D P₂) :
    (psiIstr Ψ P₁ P₂ h₁ h₂).Full where
  map_surjective {A B} g := by
    obtain ⟨f, hf⟩ := Ψ.map_surjective (g.hom : Ψ.obj A.obj ⟶ Ψ.obj B.obj)
    exact ⟨ObjectProperty.homMk f, InducedCategory.hom_ext hf⟩

/-- ★★**`Ψ^istr` は本質的全射** —— `Ψ` の本質的全射性 ＋
**isotropy が反射される**こと(第 1 主張の `mpr`)。 -/
theorem psiIstr_essSurj (h₁ : IsOfQuasiIsotropicType C P₁)
    (h₂ : IsOfQuasiIsotropicType D P₂) :
    (psiIstr Ψ P₁ P₂ h₁ h₂).EssSurj where
  mem_essImage Z := by
    obtain ⟨A₀, ⟨ε⟩⟩ : ∃ A₀ : C, Nonempty (Ψ.obj A₀ ≅ Z.obj) :=
      ⟨Ψ.objPreimage Z.obj, ⟨Ψ.objObjPreimageIso Z.obj⟩⟩
    have hiso₂ : IsIsotropic P₂ (Ψ.obj A₀) := isIsotropic_of_iso P₂ ε Z.property
    have hiso₁ : IsIsotropic P₁ A₀ :=
      (thm_3_4_i_isotropic Ψ P₁ P₂ h₁ h₂ A₀).mpr hiso₂
    exact ⟨⟨A₀, hiso₁⟩, ⟨ObjectProperty.isoMk _ ε⟩⟩

/-- ★★★★**[FrdI] Theorem 3.4, (i)** —— `Ψ^istr` は**圏同値**。

原文 (FrdI p.62):
> the horizontal arrows are equivalences of categories]. Finally, if D
-/
theorem psiIstr_isEquivalence (h₁ : IsOfQuasiIsotropicType C P₁)
    (h₂ : IsOfQuasiIsotropicType D P₂) :
    (psiIstr Ψ P₁ P₂ h₁ h₂).IsEquivalence := by
  haveI := psiIstr_faithful Ψ P₁ P₂ h₁ h₂
  haveI := psiIstr_full Ψ P₁ P₂ h₁ h₂
  haveI := psiIstr_essSurj Ψ P₁ P₂ h₁ h₂
  exact { }

/-- ★★★**`Ψ^istr` と包含関手の可換性** —— 図式の 1-可換性の**核**。

★★これは **`Iso.refl` で済む**(構成から `ι₂ ∘ Ψ^istr` と `Ψ ∘ ι₁` は
対象・射ともに同じものだから)。★原文が「1-commutative」と書く所が、
**包含関手の側では恒等的に可換**になっている。 -/
def psiIstr_comp_iota (h₁ : IsOfQuasiIsotropicType C P₁)
    (h₂ : IsOfQuasiIsotropicType D P₂) :
    psiIstr Ψ P₁ P₂ h₁ h₂ ⋙ (isotropicProp P₂).ι ≅ (isotropicProp P₁).ι ⋙ Ψ :=
  Iso.refl _

/-! ### ★★★残り —— 図式の 1-可換性と rigidity(設計)

★上の `psiIstr_comp_iota` は**包含関手**の側の可換性である。
原文の図式はもう一方、**isotropification 関手**の側:

  `Ψ^istr ∘ isotropification₁ ≅ isotropification₂ ∘ Ψ`

★★**筋(測定済み)**: どちらも**左随伴**である。
- `isotropification₁ ⊣ ι₁`(`isotropificationAdj`)、`Ψ^istr` は圏同値
- `isotropification₂ ⊣ ι₂`、`Ψ` は圏同値

★したがって `Adjunction.comp` で両者を左随伴として書き、
★★**右随伴が一致すること**を `psiIstr_comp_iota`(上)から出せば、
`Adjunction.leftAdjointUniq` で自然同型が得られる。

★**一度実装を試みて差し戻した**(2026-08-17)。詰まったのは随伴の側ではなく
**右随伴どうしの同型 `Ψ^istr⁻¹ ⋙ ι₁ ≅ ι₂ ⋙ Ψ⁻¹` の構成**である
——`psiIstr ⋙ ι₂ = ι₁ ⋙ Ψ` は**定義的に等しい**のに、
そこから逆関手側の同型を作るには単位・余単位の whisker 計算が要る。
★**次段の候補 2 つ**:
1. 上の随伴の筋を、whisker 計算を丁寧に書いて通す
2. ★**hull の一意性から直接** `NatIso.ofComponents` で作る
   —— どちらの対象も `Ψ.obj A` の hull だから(`isotropicHull_map`)。
   自然性は hull の普遍性の一意性から出る。**こちらの方が短い見込み**

★**1-一意性**も同じ随伴の一意性から出る。

★★残る `rigidity` は別筋である。原文 (FrdI p.63):
`Proposition 1.13, (ii)` を使い、**Frobenius-trivial 対象の
base-identity 自己射(任意の Frobenius 次数)に関する関手性**から
`α^d = α`(すべての `d`)を出して `α` を自明にする。
★Frobenius-normalized 型が効く箇所である。
-/

end Isotropic

end ABC3.Found.FrdI
