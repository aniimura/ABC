/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.CommonDegree
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★被覆を `Fin (N+1)` で並べ直す —— `ψ` の入力を作る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か —— 段 E2 と `§9-911` を繋ぐ配管

`§9-CommonDegree` の `exists_common_degree_cover` は

    「ある `L > 0` と**有限集合** `S ⊆ Γ(X, M^{⊗L})` があって `⨆_{u ∈ S} X_u = X`」

を与える。★一方 `§9-911` の `globalToProj` が要求するのは

    `s : Fin (N+1) → Γ(X, M^{⊗L})` と `⨆_i X_{s_i} = ⊤`

という**添字づけられた**形である。★★本ファイルがその変換を行う。

## ★機構 —— 全単射は要らない

★★**`Fin (N+1) → S` は全射でありさえすればよい**（`⨆` が増えるだけなら困らない）。
そこで `N ≔ |S|` と取り、

    `s k ≝ if k < |S| then (S の k 番目) else u₀`

と置く（`u₀` は `S` の任意の元）。★これで添字が `Fin (|S|+1)` になり、
`S` のすべての元が像に入る。
★★`|S| = 0` を避けるために `X` が**空でない**ことを使う
——`S = ∅` なら `⨆ = ∅` で被覆にならないからである。

## ★これで何が繋がったか

★★★`IsAmple M` と `X` の点 1 つから

    `∃ L > 0, ∃ N, ∃ s : Fin (N+1) → Γ(X, M^{⊗L}), ⨆_i X_{s_i} = ⊤`

が出る。★これは `globalToProj`（`§9-911`）と `isImmersion_globalToProj`（`§9-913`）の
**被覆仮定 `hcov` そのもの**である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★★★★有限被覆を `Fin (N+1)` で並べ直す -/

/-- ★★★★★★★**有限被覆を `Fin (N+1)` で並べ直す**。

★全単射は要らない——`Fin (N+1) → S` が**全射**でありさえすればよい。
`N ≔ |S|` と取り、範囲外の添字には `S` の任意の元を割り当てる。 -/
theorem exists_fin_cover (M : X.PresheafOfModules)
    (S : Set ((M.obj (op ⊤)) : Type)) (hfin : S.Finite) (hne : S.Nonempty)
    (hcov : (⨆ u ∈ S, (nonVanishing M u : Set X)) = Set.univ) :
    ∃ (N : ℕ) (s : Fin (N + 1) → ((M.obj (op ⊤)) : Type)),
      (⨆ i, nonVanishing M (s i)) = ⊤ := by
  classical
  obtain ⟨u₀, hu₀⟩ := hne
  set F := hfin.toFinset with hF
  set n := F.card with hn
  set e := (Finset.equivFin F).symm with he
  refine ⟨n, fun k =>
    if h : (k : ℕ) < n then (e ⟨(k : ℕ), h⟩ : ((M.obj (op ⊤)) : Type)) else u₀, ?_⟩
  refine eq_top_iff.2 (fun x _ => ?_)
  have hx : x ∈ (⨆ u ∈ S, (nonVanishing M u : Set X)) := by rw [hcov]; trivial
  simp only [Set.iSup_eq_iUnion, Set.mem_iUnion, exists_prop] at hx
  obtain ⟨u, huS, hxu⟩ := hx
  have huF : u ∈ F := hfin.mem_toFinset.2 huS
  set k := (Finset.equivFin F) ⟨u, huF⟩ with hk
  have hklt : (k : ℕ) < n := k.2
  refine TopologicalSpace.Opens.mem_iSup.2 ⟨⟨(k : ℕ), by omega⟩, ?_⟩
  simp only [hklt, dif_pos]
  have hval : e ⟨(k : ℕ), hklt⟩ = (⟨u, huF⟩ : F) := by
    rw [he]; simp [hk]
  rw [hval]
  exact hxu

/-! ## ★★★★★★★★★★`IsAmple` から直接 -/

/-- ★★★★★★★★★★**`IsAmple` から `ψ` の入力を作る**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `IsAmple M` ⟹ `∃ L > 0, ∃ N, ∃ s : Fin (N+1) → Γ(X, M^{⊗L}), ⨆_i X_{s_i} = ⊤`

★これが `globalToProj`（`§9-911`）と `isImmersion_globalToProj`（`§9-913`）の
**被覆仮定 `hcov` そのもの**である。
★★`X` の点 1 つを受けるのは `S = ∅` を排除するためである。 -/
theorem exists_fin_cover_of_isAmple (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (h : IsAmple M) (x₀ : (X : Type)) :
    ∃ L : ℕ, 0 < L ∧ ∃ (N : ℕ)
      (s : Fin (N + 1) → ((presheafTensorPow M L).obj (op ⊤) : Type)),
      (⨆ i, nonVanishing (presheafTensorPow M L) (s i)) = ⊤ := by
  obtain ⟨L, hL, S, hSfin, hScov⟩ := exists_common_degree_cover M hM h
  refine ⟨L, hL, ?_⟩
  have hne : S.Nonempty := by
    by_contra hemp
    rw [Set.not_nonempty_iff_eq_empty] at hemp
    subst hemp
    simp only [Set.mem_empty_iff_false, iSup_false, iSup_bot] at hScov
    have hcontra : x₀ ∈ (⊥ : Set X) := by rw [hScov]; trivial
    simp at hcontra
  exact exists_fin_cover _ S hSfin hne hScov

/-! ## ★出典の紐付け(`.src`) -/

def exists_fin_cover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(有限被覆を Fin (N+1) で並べ直す)",
    sectionId := "genell-prop-1-4" }

def exists_fin_cover_of_isAmple.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(IsAmple から ψ の入力を作る)",
    sectionId := "genell-prop-1-4" }

def exists_fin_cover_of_isAmple.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_common_degree_cover(ample なら共通次数の有限個の切断で覆える)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_common_degree_cover") 2,
    .implicitStep
      ("★全単射は要らない——Fin (N+1) → S が**全射**でありさえすればよい。" ++
       "N ≔ |S| と取り、範囲外の添字には S の任意の元を割り当てる") 1,
    .implicitStep
      ("★★これで globalToProj(§9-911)と isImmersion_globalToProj(§9-913)の" ++
       "被覆仮定 hcov が IsAmple から出る。" ++
       "★残るのは haff(チャートがアフィン)と hsurj(チャート写像が全射)であり、" ++
       "どちらも段 E3 の内容である") 4 ]

end ABC3.Found.GenEll
