import ABC3.Found.FrdI.Def31

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

/-! ### ★★★段 5 の設計(測定済み・実装待ち)

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

end ABC3.Found.FrdI
