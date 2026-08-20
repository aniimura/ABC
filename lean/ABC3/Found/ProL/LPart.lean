/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Found.ProL.Defs
import ABC3.Found.ProL.ProfiniteLimit
import Mathlib.Topology.Algebra.Group.ClosedSubgroup
import Mathlib.Topology.Algebra.ClopenNhdofOne

/-!
# `M[l]` —— 可換副有限群の pro-`l` 部分(チェーン `prol` の葉 `M-l-part`)

最終目標は [FrdI] `Definition 2.8, (ii)`(物理 p.52)。

原文 (FrdI p.52):
> Thus, M decomposes as a direct product of pro-l groups M [l], where l varies over

## ★★測って選んだ構成(2026-08-19)

台帳は `M[l] := lim_U (M/U)[l]` と書いていたが、
★**極限として組むより `M` の閉部分群として取るほうが軽い** ——
`ProfiniteGrp.ofClosedSubgroup` が在庫にあるからである。

`M[l] := {x : M | すべての有限商 M/U で x の像が l-準素}`

★★**可換性は型クラスで持つ**(`[CommGroup M]`)。仮定 `hcomm : ∀ a b, a*b = b*a` の形だと
下流の `primProj`(`CommGroup` を要求する)を使うたびに `letI` が要り、
**同じ instance を作り直す羽目になる**からである。
`ProfiniteGrp` の側へは最後に一度だけ橋を架ける。
-/

namespace ABC3.Found.ProL

open CategoryTheory

universe u

section CommProfinite

variable {M : Type u} [CommGroup M] [TopologicalSpace M] [IsTopologicalGroup M]
  [CompactSpace M] [TotallyDisconnectedSpace M]

variable (M) in
/-- ★`M[l]` の台 —— すべての有限商で `l`-準素な元。 -/
def lPartCarrier (l : ℕ) : Set M :=
  {x : M | ∀ U : OpenNormalSubgroup M, ∃ k : ℕ,
    (QuotientGroup.mk x : M ⧸ U.toSubgroup) ^ l ^ k = 1}

theorem mem_lPartCarrier {l : ℕ} {x : M} :
    x ∈ lPartCarrier M l ↔ ∀ U : OpenNormalSubgroup M, ∃ k : ℕ,
      (QuotientGroup.mk x : M ⧸ U.toSubgroup) ^ l ^ k = 1 := Iff.rfl

variable (M) in
/-- ★★**`M[l]`** —— `M` の部分群。 -/
def lPart (l : ℕ) : Subgroup M where
  carrier := lPartCarrier M l
  one_mem' := fun U => ⟨0, by simp⟩
  mul_mem' {x y} hx hy := by
    intro U
    obtain ⟨k, hk⟩ := hx U
    obtain ⟨m, hm⟩ := hy U
    refine ⟨k + m, ?_⟩
    have hx' : (QuotientGroup.mk x : M ⧸ U.toSubgroup) ^ l ^ (k + m) = 1 := by
      rw [pow_add, pow_mul, hk, one_pow]
    have hy' : (QuotientGroup.mk y : M ⧸ U.toSubgroup) ^ l ^ (k + m) = 1 := by
      rw [pow_add, pow_mul', hm, one_pow]
    rw [QuotientGroup.mk_mul, mul_pow, hx', hy', one_mul]
  inv_mem' {x} hx := by
    intro U
    obtain ⟨k, hk⟩ := hx U
    exact ⟨k, by rw [QuotientGroup.mk_inv, inv_pow, hk, inv_one]⟩

variable (M) in
/-- ★★`M[l]` は閉集合 —— 各有限商での条件は離散空間への連続写像の逆像だから。 -/
theorem lPartCarrier_isClosed (l : ℕ) : IsClosed (lPartCarrier M l) := by
  have hset : lPartCarrier M l = ⋂ U : OpenNormalSubgroup M,
      (fun x : M => (QuotientGroup.mk x : M ⧸ U.toSubgroup))
        ⁻¹' {z : M ⧸ U.toSubgroup | ∃ k : ℕ, z ^ l ^ k = 1} := by
    ext x
    simp only [lPartCarrier, Set.mem_setOf_eq, Set.mem_iInter, Set.mem_preimage]
  rw [hset]
  refine isClosed_iInter (fun U => ?_)
  haveI : DiscreteTopology (M ⧸ U.toSubgroup) :=
    inferInstanceAs (DiscreteTopology (M ⧸ U.toOpenSubgroup.toSubgroup))
  exact (isClosed_discrete _).preimage continuous_quotient_mk'

instance lPart_compactSpace (l : ℕ) : CompactSpace (lPart M l) :=
  isCompact_iff_compactSpace.mp ((lPartCarrier_isClosed M l).isCompact)

variable (M) in
/-- ★★★**`M[l]` は副有限群**。 -/
def lPartGrp (l : ℕ) : ProfiniteGrp.{u} :=
  ProfiniteGrp.of (lPart M l)

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`M[l]` は pro-`l` 群**。

★`M[l]` の開正規部分群 `V` に対し、`M` の開正規部分群 `U` で
`M[l] ∩ U ⊆ V` となるものを取る(`exist_openNormalSubgroup_sub_open_nhds_of_one`)。
★あとは `M/U` での `l`-準素性をそのまま運ぶ。 -/
theorem isProL_lPartGrp (l : ℕ) : IsProL l (lPartGrp M l) := by
  intro V y
  obtain ⟨x, rfl⟩ := QuotientGroup.mk_surjective y
  set xs : (lPart M l) := x with hxs
  obtain ⟨W, hWopen, hWpre⟩ : ∃ W : Set M, IsOpen W ∧
      Subtype.val ⁻¹' W = ((V : Subgroup (lPartGrp M l)) : Set (lPartGrp M l)) := by
    obtain ⟨W, hWopen, hWeq⟩ := isOpen_induced_iff.mp V.isOpen'
    exact ⟨W, hWopen, hWeq⟩
  have h1W : (1 : M) ∈ W := by
    have h0 : (1 : lPartGrp M l) ∈ ((V : Subgroup (lPartGrp M l)) : Set (lPartGrp M l)) :=
      V.toSubgroup.one_mem
    rw [← hWpre] at h0
    exact h0
  obtain ⟨U, hU⟩ := ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one
    (G := M) hWopen h1W
  obtain ⟨k, hk⟩ := (mem_lPartCarrier (l := l) (M := M)).mp xs.2 U
  refine ⟨k, ?_⟩
  rw [← QuotientGroup.mk_pow, QuotientGroup.eq_one_iff]
  have hxU : (xs.1 : M) ^ l ^ k ∈ U.toSubgroup := by
    rw [← QuotientGroup.mk_pow, QuotientGroup.eq_one_iff] at hk
    exact hk
  have hxW : (xs.1 : M) ^ l ^ k ∈ W := hU hxU
  have h2 : (xs ^ l ^ k) ∈ Subtype.val ⁻¹' W := hxW
  rw [hWpre] at h2
  exact h2

def lPartGrp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52, item := "Definition 2.8, (ii) — M[l] の構成",
    sectionId := "frdi-def-2-8" }

end CommProfinite

end ABC3.Found.ProL
