/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Topology.Algebra.ClopenNhdofOne
import Mathlib.Topology.Algebra.Group.Basic
import Mathlib.GroupTheory.Index
import Mathlib.GroupTheory.ResiduallyFinite
import Mathlib.Tactic.Group
import ABC3.Meta.Claim

/-!
# 副有限群は residually finite —— `Theorem 6.2, (iv)` の核(`Found`、`sorry` 無し)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.111。

原文 (FrdI p.111):
> nioid C is of isotropic, standard, and birationally Frobenius-normalized

## ★★なにを閉じたか

原文 `Theorem 6.2, (iv)` は `Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)` から
`C` の **Frobenius-slim / slim / Div-slim** を「(formally)」の 1 語で出し、
[Mzk7] `Corollary 1.1.6` へ送る。その **1 語の中身**が

> 副有限群 `G` は residually finite である
> ——`g ≠ 1` なら `g ∉ N` となる**開正規**部分群 `N` が取れる

であり、これを `sorry` 無しで閉じたのが本ファイルである。

## ★★★測定 —— mathlib に「開正規部分群」版は無かった(2026-08-25)

mathlib が持っているのは

* `IsTopologicalGroup.exist_openSubgroup_sub_clopen_nhds_of_one`
  …… clopen 近傍の中に**開部分群**を取る(**正規とは限らない**)
* `Subgroup.normalCore` / `normalCore_le` / `normalCore_normal`
  …… 正規核(**開性は言わない**)

★両者を繋ぐ「**開部分群の正規核は開**」が無かったので、それを
`isOpen_normalCore_of_openSubgroup` として書く。中身は

```
H.normalCore = ⋂_{q : G ⧸ H} (a ↦ (out q)⁻¹ * a * out q)⁻¹ (H)
```

——★**添字が有限**であることが要点で、それは `G` がコンパクトで `H` が開
(したがって `G ⧸ H` が離散コンパクト＝有限)から出る。

★★共役 `b H b⁻¹` は左剰余類 `bH` だけで決まる、と書きたくなるが、
`normalCore = {a | ∀ b, b a b⁻¹ ∈ H}` を素直に添字付けると**右**剰余類になる。
そこで `b ↦ b⁻¹` で読み替えて `{a | ∀ b, b⁻¹ a b ∈ H}` の形にすると、
`(bh)⁻¹ a (bh) = h⁻¹ (b⁻¹ a b) h` なので `h ∈ H` の分が `H` の中で吸収され、
**左**剰余類 `G ⧸ H` で添字付けできる。これが下の `hset` の中身である。
-/

namespace ABC3.Found.FrdI.Profinite

open ABC3.Meta

/-! ## ★1. 開部分群の正規核は開 -/

/-- ★★★★★**開部分群の正規核は開**(コンパクト位相群)。

`H.normalCore = ⋂_{q : G ⧸ H} (a ↦ (out q)⁻¹ a (out q))⁻¹ (H)` と書き直す。
`G` コンパクト・`H` 開なので `G ⧸ H` は有限、各項は連続写像による開集合の逆像。 -/
theorem isOpen_normalCore_of_openSubgroup {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] (H : OpenSubgroup G) :
    IsOpen (((H : Subgroup G).normalCore : Subgroup G) : Set G) := by
  haveI : (H : Subgroup G).FiniteIndex := Subgroup.finiteIndex_of_finite_quotient
  have hset : (((H : Subgroup G).normalCore : Subgroup G) : Set G)
      = ⋂ (q : G ⧸ (H : Subgroup G)),
        (fun a => (Quotient.out q)⁻¹ * a * (Quotient.out q)) ⁻¹' ((H : Subgroup G) : Set G) := by
    ext a
    simp only [Set.mem_iInter, Set.mem_preimage, SetLike.mem_coe]
    constructor
    · intro ha q
      have := ha (Quotient.out q)⁻¹
      rwa [inv_inv] at this
    · intro ha b
      -- ★`out ⟦b⁻¹⟧` は `b⁻¹` と同じ左剰余類にいるので `b⁻¹ * h` の形
      obtain ⟨h, hh, hout⟩ : ∃ h ∈ (H : Subgroup G),
          Quotient.out (QuotientGroup.mk b⁻¹ : G ⧸ (H : Subgroup G)) = b⁻¹ * h := by
        refine ⟨b * Quotient.out (QuotientGroup.mk b⁻¹ : G ⧸ (H : Subgroup G)), ?_, by group⟩
        have h1 : (QuotientGroup.mk (Quotient.out (QuotientGroup.mk b⁻¹ : G ⧸ (H : Subgroup G)))
            : G ⧸ (H : Subgroup G)) = QuotientGroup.mk b⁻¹ := Quotient.out_eq _
        rw [QuotientGroup.eq] at h1
        simpa using (H : Subgroup G).inv_mem h1
      have hq := ha (QuotientGroup.mk b⁻¹)
      rw [hout] at hq
      -- ★`h` の分は `H` の中で吸収される
      have hc : h⁻¹ * (b * a * b⁻¹) * h ∈ (H : Subgroup G) := by
        convert hq using 1; group
      have h2 := (H : Subgroup G).mul_mem
        ((H : Subgroup G).mul_mem hh hc) ((H : Subgroup G).inv_mem hh)
      convert h2 using 1; group
  rw [hset]
  refine isOpen_iInter_of_finite fun q => ?_
  exact IsOpen.preimage ((continuous_const.mul continuous_id).mul continuous_const) H.isOpen'

/-! ## ★2. 副有限群は residually finite -/

/-- ★★★★★★**副有限群は residually finite** ——
`1` でない元は、ある**開正規**部分群の外にある。

★これが `Theorem 6.2, (iv)` の「(formally)」の中身である。

★仮定「コンパクト・Hausdorff・完全不連結な位相群」が副有限性の言い換え。
`isTopologicalBasis_isClopen` で `1 ∈ W ⊆ {g}ᶜ` なる clopen を取り、
`exist_openSubgroup_sub_clopen_nhds_of_one` でその中の開部分群 `H` を取り、
`H.normalCore` を返す。 -/
theorem exists_open_normal_notMem {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] [T2Space G] [TotallyDisconnectedSpace G]
    (g : G) (hg : g ≠ 1) :
    ∃ N : Subgroup G, N.Normal ∧ IsOpen (N : Set G) ∧ g ∉ N := by
  obtain ⟨W, hWclopen, hWmem, hWsub⟩ :=
    isTopologicalBasis_isClopen.exists_subset_of_mem_open (a := (1 : G)) (u := {g}ᶜ)
      (by simpa using hg.symm) isOpen_compl_singleton
  obtain ⟨H, hHW⟩ := IsTopologicalGroup.exist_openSubgroup_sub_clopen_nhds_of_one hWclopen hWmem
  refine ⟨(H : Subgroup G).normalCore, Subgroup.normalCore_normal _,
    isOpen_normalCore_of_openSubgroup H, ?_⟩
  intro hmem
  exact (hWsub (hHW (Subgroup.normalCore_le _ hmem))) rfl

/-- ★★**開正規部分群の交わりは自明** —— residually finite の言い換え。 -/
theorem iInf_open_normal_eq_bot {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] [T2Space G] [TotallyDisconnectedSpace G]
    (g : G) (hg : ∀ N : Subgroup G, N.Normal → IsOpen (N : Set G) → g ∈ N) :
    g = 1 := by
  by_contra hne
  obtain ⟨N, hNn, hNo, hNg⟩ := exists_open_normal_notMem g hne
  exact hNg (hg N hNn hNo)

/-! ## ★3. mathlib の `Group.ResiduallyFinite` に載せる -/

/-- ★★**開部分群は有限指数**(コンパクト位相群)。 -/
theorem finiteIndex_of_isOpen {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] (N : Subgroup G) (hNo : IsOpen (N : Set G)) :
    N.FiniteIndex := by
  haveI : Finite (G ⧸ N) :=
    Subgroup.instFiniteQuotientOfSeparatelyContinuousMulOfCompactSpace (⟨N, hNo⟩ : OpenSubgroup G)
  exact Subgroup.finiteIndex_of_finite_quotient

/-- ★★★★★★**副有限群は `Group.ResiduallyFinite`** ——
mathlib の型クラスに載せる。これで `Remark 3.1.2`
(`isFrobeniusSlim_of_residuallyFinite`)が直接使える。

★mathlib には「コンパクト Hausdorff 完全不連結位相群 ⟹ `ResiduallyFinite`」の
インスタンスが**無い**(2026-08-25 実測)。 -/
theorem residuallyFinite_of_profinite {G : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [CompactSpace G] [T2Space G] [TotallyDisconnectedSpace G] :
    Group.ResiduallyFinite G := by
  rw [Group.residuallyFinite_iff_exists_finiteIndex]
  intro g hg
  obtain ⟨N, _, hNo, hNg⟩ := exists_open_normal_notMem g hg
  exact ⟨N, finiteIndex_of_isOpen N hNo, hNg⟩

/-- ★★★**単射準同型に沿って `ResiduallyFinite` は下りる**。

★mathlib には部分群版(`Group.instResiduallyFiniteSubtypeMemSubgroup`)しか無い
(2026-08-25 実測)ので、`Aut(𝒟_A → 𝒟) ↪ G` の形に当てるためにこちらを書く。 -/
theorem residuallyFinite_of_injective {G H : Type*} [Group G] [Group H]
    [Group.ResiduallyFinite G] (f : H →* G) (hf : Function.Injective f) :
    Group.ResiduallyFinite H := by
  refine Group.residuallyFinite_of_forall_exists_finite_monoidHom (fun h hh => ?_)
  have hfh : f h ≠ 1 := fun hc => hh (hf (by simpa using hc))
  obtain ⟨N, hNg⟩ := Group.exists_finiteIndexNormalSubgroup_notMem (f h) hfh
  haveI : N.toSubgroup.Normal := N.isNormal'
  haveI : N.toSubgroup.FiniteIndex := N.isFiniteIndex'
  haveI : Finite (G ⧸ N.toSubgroup) := N.toSubgroup.finite_quotient_of_finiteIndex
  refine ⟨G ⧸ N.toSubgroup, inferInstance, inferInstance,
    (QuotientGroup.mk' N.toSubgroup).comp f, ?_⟩
  have hmem : f h ∉ N.toSubgroup := hNg
  simpa [QuotientGroup.eq_one_iff] using hmem

/-- ★★**同型に沿って `ResiduallyFinite` は移る**。 -/
theorem residuallyFinite_of_mulEquiv {G H : Type*} [Group G] [Group H]
    [Group.ResiduallyFinite G] (e : H ≃* G) : Group.ResiduallyFinite H :=
  residuallyFinite_of_injective (e : H →* G) e.injective

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_open_normal_notMem.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — 副有限群は residually finite",
    sectionId := "frdi-thm-6-2" }

def exists_open_normal_notMem.needs : List ProofObligation :=
  [ .citation "[mathlib]" "IsTopologicalGroup.exist_openSubgroup_sub_clopen_nhds_of_one"
      (.inMathlib "IsTopologicalGroup.exist_openSubgroup_sub_clopen_nhds_of_one") 111,
    .citation "[mathlib]" "isTopologicalBasis_isClopen(完全不連結コンパクト Hausdorff)"
      (.inMathlib "isTopologicalBasis_isClopen") 111,
    .citation "[mathlib]" "Subgroup.normalCore / normalCore_normal / normalCore_le"
      (.inMathlib "Subgroup.normalCore") 111,
    .citation "[mathlib]" "開部分群の正規核が開であること"
      (.absent "Subgroup.normalCore 周辺と OpenSubgroup 周辺の宣言を列挙(2026-08-25)。normalCore の開性・閉性を述べる宣言は無い") 111,
    .derivation
      "H.normalCore = ⋂_{q : G ⧸ H} (a ↦ (out q)⁻¹ a (out q))⁻¹ (H) と書き直す。G ⧸ H が有限なので有限交叉" 111,
    .implicitStep
      "★原文は「(formally)」の 1 語で畳み、[Mzk7] Corollary 1.1.6 へ送る" 111 ]

def isOpen_normalCore_of_openSubgroup.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — 開部分群の正規核は開",
    sectionId := "frdi-thm-6-2" }

def isOpen_normalCore_of_openSubgroup.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Subgroup.finiteIndex_of_finite_quotient"
      (.inMathlib "Subgroup.finiteIndex_of_finite_quotient") 111,
    .derivation "共役 b⁻¹ H b は左剰余類 bH だけで決まる((bh)⁻¹ a (bh) = h⁻¹ (b⁻¹ a b) h)" 111 ]

def iInf_open_normal_eq_bot.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — 開正規部分群の交わりは自明",
    sectionId := "frdi-thm-6-2" }

def iInf_open_normal_eq_bot.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_open_normal_notMem"
      (.inProject "ABC3" "ABC3.Found.FrdI.Profinite.exists_open_normal_notMem") 111 ]

def residuallyFinite_of_profinite.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — Since G is profinite, hence, in particular, residually finite",
    sectionId := "frdi-thm-6-2" }

def residuallyFinite_of_profinite.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Group.ResiduallyFinite / residuallyFinite_iff_exists_finiteIndex"
      (.inMathlib "Group.ResiduallyFinite") 111,
    .citation "[mathlib]" "副有限群が ResiduallyFinite であるインスタンス"
      (.absent "Group.ResiduallyFinite のインスタンスは Finite / Subgroup / Prod の 3 つのみ(2026-08-25 実測)。位相からのものは無い") 111,
    .citation "[ABC3]" "exists_open_normal_notMem"
      (.inProject "ABC3" "ABC3.Found.FrdI.Profinite.exists_open_normal_notMem") 111 ]

def residuallyFinite_of_injective.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iv) — it follows formally that Z, Z_G(H) are also residually finite",
    sectionId := "frdi-thm-6-2" }

def residuallyFinite_of_injective.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Group.residuallyFinite_of_forall_exists_finite_monoidHom"
      (.inMathlib "Group.residuallyFinite_of_forall_exists_finite_monoidHom") 111,
    .derivation "単射準同型で押し出し、G の有限商へ合成する" 111 ]

end ABC3.Found.FrdI.Profinite
