import ABC3.Found.PGC.LubinTateReciprocitySurjective
import ABC3.Found.PGC.LubinTateActionEquivariance

/-!
# `K.carrier⟮x⟯` は `K.carrier` 上正規(`sorry` 無し)

節目(5)(射影極限 `Gal(L_π/K)≅𝒪_K^×`)の**さらに先**——古典的Lubin-Tate
理論が実際に主張する`Gal(K_π/K)≅𝒪_K^×`(`K_π:=K(Λ_∞)`、完全な代数
閉包ではない)へ向けた布石。`K_π`を`⨆n,K.carrier⟮(psiGenSeq n).pt⟯`
として構成する際、`IntermediateField.normal_iSup`(mathlib)で
`Normal K.carrier K_π`を得るには**各段**`K.carrier⟮x⟯`(`x`は
`ψ_n`の根)が`Normal K.carrier`であることが要る——これを示す。

## 証明の構造

`x`の最小多項式は`ψ_n`自身(`minpoly.eq_of_irreducible_of_monic`、
`Found/PGC/LubinTateReciprocitySurjective.lean`と同じ手筋)。`ψ_n`の
根はすべて`K.carrier⟮x⟯`に収まる(`iteratedLubinTatePsiTorsionPoints_
subset_adjoin`、既出)ので、`K.carrier⟮x⟯`は`x`の最小多項式の
**分裂体**(`IntermediateField.adjoin_rootSet_isSplittingField`で
`adjoin(rootSet)`として作った分裂体が、根がすべて`K.carrier⟮x⟯`に
収まる+`x`自身が根であることから`K.carrier⟮x⟯`と一致する)。
`Normal.of_isSplittingField`で結論する。

★別ルートとの関係: `Found/PGC/LubinTateReciprocityIsomorphism.lean::
galoisReciprocityEquiv`(既出、`Gal(K.carrier⟮x⟯/K.carrier)≃*
(𝒪_K/π^n)^×`)からも、`#Aut(K(x)/K)=[K(x):K]`(自己同型の個数が
拡大次数に一致)経由でNormal性を導ける見込みだが、そちらは
`IsGalois.of_card_aut_eq_finrank`のような追加の橋渡しが要る。
ここでは分裂体からの直接証明を採った——`K_π`のsup-normal性の構成
にすぐ使える形。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- **`K.carrier⟮x⟯`は`K.carrier`上正規**——`x`が`ψ_n`の根なら、
`K.carrier⟮x⟯`は`x`の最小多項式(`=ψ_n`)の分裂体に一致する。 -/
theorem normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  haveI := valuationRing_isDVR K
  haveI := IsAlgClosure.normal K.carrier K.closure
  set psiK : Polynomial K.carrier :=
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map (algebraMap (𝒪[K.carrier]) K.carrier)
    with hpsiK
  have hmap : algebraMap (𝒪[K.carrier]) K.closure =
      (algebraMap K.carrier K.closure).comp (algebraMap (𝒪[K.carrier]) K.carrier) :=
    IsScalarTower.algebraMap_eq (𝒪[K.carrier]) K.carrier K.closure
  have hrw : (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) =
      psiK.map (algebraMap K.carrier K.closure) := by
    rw [hmap, ← Polynomial.map_map, hpsiK]
  -- `iteratedLubinTatePsiTorsionPoints`の定義(`(mapされたψ_n).roots.toFinset`)を
  -- mathlibの`Polynomial.rootSet`(`Polynomial.mem_rootSet'`)の言葉へ翻訳する。
  have hiff : ∀ z : K.closure, z ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn ↔
      z ∈ psiK.rootSet K.closure := by
    intro z
    rw [iteratedLubinTatePsiTorsionPoints, Multiset.mem_toFinset, hrw, Polynomial.mem_roots',
      Polynomial.mem_rootSet']
    constructor
    · intro ⟨h1, h2⟩
      refine ⟨h1, ?_⟩
      rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] at h2
      exact h2
    · intro ⟨h1, h2⟩
      refine ⟨h1, ?_⟩
      rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def]
      exact h2
  have hsplits : (psiK.map (algebraMap K.carrier K.closure)).Splits :=
    IsAlgClosed.splits (Polynomial.map (algebraMap K.carrier K.closure) psiK)
  have hsplitfield := IntermediateField.adjoin_rootSet_isSplittingField hsplits
  -- `adjoin(rootSet)`は、根がすべて`K.carrier⟮x⟯`に収まる(⊆)ことと
  -- `x`自身が根であること(⊇)から、`K.carrier⟮x⟯`そのものに一致する。
  have heq : IntermediateField.adjoin K.carrier (psiK.rootSet K.closure) =
      IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
    apply le_antisymm
    · rw [IntermediateField.adjoin_le_iff]
      intro z hz
      exact iteratedLubinTatePsiTorsionPoints_subset_adjoin K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
        z ((hiff z).mpr hz)
    · rw [IntermediateField.adjoin_le_iff]
      intro z hz
      simp only [Set.mem_singleton_iff] at hz
      rw [hz]
      exact IntermediateField.subset_adjoin K.carrier (psiK.rootSet K.closure) ((hiff x).mp hxψ)
  rw [heq] at hsplitfield
  exact Normal.of_isSplittingField psiK

end ABC3.Found.PGC
