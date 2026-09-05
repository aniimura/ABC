import Mathlib.RingTheory.Henselian
import Mathlib.FieldTheory.Separable
import Mathlib.Algebra.Polynomial.Splits

/-!
# Henselian 局所環: 剰余体で分裂する分離モニック多項式は、商体でも分裂する

論文にも我々のモデルにも固有でない、**一般の**結果。mathlib へ出せる形で書く
(`Found/ResidueFieldFinite.lean`・`Found/FiniteFieldIrreducible.lean` と同じ位置づけ)。

## 何を言っているか

`A` を Henselian 局所整域、`F : A[X]` をモニックとする。`F` の剰余体への像 `F̄` が

* **分離的**(重根を持たない)で、
* 剰余体 `𝓀_A` の中で**分裂する**(1 次式の積になる)

ならば、`F` 自身も `A` を埋め込む任意の体 `L` の中で**分裂する**。

## なぜ要るか

`Found/PGC/UnramifiedExtension.lean` で、不分岐拡大 `K(x)/K` が **normal**
(したがって標数 0 なので Galois)であることを示すのに使う——不分岐拡大の
定義多項式 `f` は、剰余体では既約分離多項式 `ḡ` になり、`𝓀_{K(x)}` は
`𝓀` の有限次拡大なので **normal**、よって `ḡ` は `𝓀_{K(x)}` で分裂する。
その分裂を Hensel で `𝒪_{K(x)}` へ持ち上げれば `f` が `K(x)` で分裂する。

## 証明の骨格

1. `F̄` の相異なる根はちょうど `deg F` 個(分離+分裂)。
2. 各根 `b` に対し、`F̄'(b) ≠ 0`(分離性)なので Hensel が
   `F(a)=0`・`residue a = b` なる `a ∈ A` を与える
   (`exists_root_of_residue_root`)。
3. すなわち `residue` は `F.roots` から `F̄.roots` への**全射**。
   よって `F` は `A` の中に `deg F` 個以上の相異なる根を持つ。
4. `ι` は単射なので `F.map ι` も `deg F` 個以上の相異なる根を持ち、
   一方 `roots.card ≤ natDegree` だから等号——`splits_iff_card_roots`。
-/

namespace ABC3.Found

open Polynomial

/-- **Henselian 局所環における単根の持ち上げ**——剰余体の単根
(`F̄(b) = 0` かつ `F̄'(b) ≠ 0`)は `A` の根に持ち上がり、その剰余は `b`。 -/
theorem exists_root_of_residue_root {A : Type*} [CommRing A] [HenselianLocalRing A]
    (F : Polynomial A) (hF : F.Monic) (b : IsLocalRing.ResidueField A)
    (hb : Polynomial.eval b (F.map (IsLocalRing.residue A)) = 0)
    (hb' : Polynomial.eval b (Polynomial.derivative (F.map (IsLocalRing.residue A))) ≠ 0) :
    ∃ a : A, F.IsRoot a ∧ IsLocalRing.residue A a = b := by
  obtain ⟨a₀, ha₀⟩ := IsLocalRing.residue_surjective b
  have h1 : Polynomial.eval a₀ F ∈ IsLocalRing.maximalIdeal A := by
    rw [← IsLocalRing.residue_eq_zero_iff]
    have := Polynomial.hom_eval₂ F (RingHom.id A) (IsLocalRing.residue A) a₀
    simp only [Polynomial.eval₂_id, RingHom.comp_id] at this
    rw [this, ← Polynomial.eval_map, ha₀]
    exact hb
  have h2 : IsUnit (Polynomial.eval a₀ (Polynomial.derivative F)) := by
    refine (IsLocalRing.residue_ne_zero_iff_isUnit _).mp ?_
    have := Polynomial.hom_eval₂ (Polynomial.derivative F) (RingHom.id A)
      (IsLocalRing.residue A) a₀
    simp only [Polynomial.eval₂_id, RingHom.comp_id] at this
    rw [this, ← Polynomial.eval_map, ha₀, ← Polynomial.derivative_map]
    exact hb'
  obtain ⟨a, ha, hsub⟩ := HenselianLocalRing.is_henselian F hF a₀ h1 h2
  refine ⟨a, ha, ?_⟩
  have : IsLocalRing.residue A (a - a₀) = 0 := (IsLocalRing.residue_eq_zero_iff _).mpr hsub
  rw [map_sub, sub_eq_zero] at this
  rw [this, ha₀]

open scoped Classical in
/-- **剰余体で分裂する分離モニック多項式は、`A` を埋め込む体でも分裂する**。 -/
theorem splits_map_of_residue_splits {A : Type*} [CommRing A] [IsDomain A] [HenselianLocalRing A]
    {L : Type*} [Field L] (ι : A →+* L) (hι : Function.Injective ι)
    (F : Polynomial A) (hFm : F.Monic)
    (hsep : (F.map (IsLocalRing.residue A)).Separable)
    (hspl : (F.map (IsLocalRing.residue A)).Splits) :
    (F.map ι).Splits := by
  have hFne : F ≠ 0 := hFm.ne_zero
  have hFLne : F.map ι ≠ 0 := (hFm.map ι).ne_zero
  have hFbne : F.map (IsLocalRing.residue A) ≠ 0 := (hFm.map _).ne_zero
  have hcardb : (F.map (IsLocalRing.residue A)).roots.toFinset.card = F.natDegree := by
    rw [Multiset.card_toFinset, Multiset.dedup_eq_self.mpr (Polynomial.nodup_roots hsep)]
    exact (Polynomial.splits_iff_card_roots.mp hspl).trans (hFm.natDegree_map _)
  have hsurj : Set.SurjOn (IsLocalRing.residue A) (↑F.roots.toFinset)
      (↑(F.map (IsLocalRing.residue A)).roots.toFinset) := by
    intro b hb
    simp only [Multiset.mem_toFinset, Finset.mem_coe] at hb
    have hbroot : (F.map (IsLocalRing.residue A)).IsRoot b := (Polynomial.mem_roots hFbne).mp hb
    have hbder : Polynomial.eval b
        (Polynomial.derivative (F.map (IsLocalRing.residue A))) ≠ 0 := by
      obtain ⟨u, v, huv⟩ := hsep
      intro hcon
      have hev := congrArg (Polynomial.eval b) huv
      simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_one, hcon,
        Polynomial.IsRoot.def.mp hbroot, mul_zero, zero_add] at hev
      exact one_ne_zero hev.symm
    obtain ⟨a, ha, hres⟩ := exists_root_of_residue_root F hFm b
      (Polynomial.IsRoot.def.mp hbroot) hbder
    exact ⟨a, by simpa using (Polynomial.mem_roots hFne).mpr ha, hres⟩
  have hle1 : (F.map (IsLocalRing.residue A)).roots.toFinset.card ≤ F.roots.toFinset.card :=
    Finset.card_le_card_of_surjOn _ hsurj
  have hle2 : F.roots.toFinset.card ≤ (F.map ι).roots.toFinset.card := by
    refine Finset.card_le_card_of_injOn ι (fun a ha => ?_) (fun a _ b _ h => hι h)
    have hroot := (Polynomial.mem_roots hFne).mp (Multiset.mem_toFinset.mp ha)
    refine Multiset.mem_toFinset.mpr ((Polynomial.mem_roots hFLne).mpr ?_)
    show Polynomial.eval (ι a) (F.map ι) = 0
    rw [Polynomial.eval_map, Polynomial.eval₂_at_apply, Polynomial.IsRoot.def.mp hroot, map_zero]
  have hle3 : (F.map ι).roots.toFinset.card ≤ Multiset.card (F.map ι).roots :=
    Multiset.toFinset_card_le _
  have hle4 : Multiset.card (F.map ι).roots ≤ (F.map ι).natDegree := Polynomial.card_roots' _
  have hdeg : (F.map ι).natDegree = F.natDegree := hFm.natDegree_map _
  rw [Polynomial.splits_iff_card_roots]
  omega

end ABC3.Found
