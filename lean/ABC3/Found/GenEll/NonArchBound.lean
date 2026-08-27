/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.CompactBound
import Mathlib.NumberTheory.Padics.PadicNumbers
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure

/-!
# [GenEll] Example 1.3, (ii) —— **非アルキメデス側**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–p.6。

原文 (GenEll p.6):
> for the subset of points x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, such that for each

## ★★何が残っていたか

`CompactBound.lean` は `Example 1.3, (ii)` の**アルキメデス側**を取り、
`.src` にも「**アルキメデス側のみ——非アルキメデス側は未着手**」と書いていた。
★本ファイルがその非アルキメデス側を作る。

## ★★★原文の定義

> For each `v ∈ V^arc`（respectively, `v ∈ V^non`), let
> `K_v ⊆ X^arc`（respectively, `K_v ⊆ X(ℚ̄_v)`)
> be a nonempty `ι_X`-stable compact domain（respectively, a nonempty
> `Gal(ℚ̄_v/ℚ_v)`-stable subset whose intersection with each `X(K) ⊆ X(ℚ̄_v)`,
> for `K ⊆ ℚ̄_v` a finite extension of `ℚ_v`, is a compact domain in `X(K)`）

★★アルキメデス側と**同じ形**である —— 「`x` が定める点が `K_v` に含まれる」。
違いは点の住む場所が `X^arc` から `X(ℚ̄_p)` に替わることだけである。

## ★★★★「compact domain」を定義に入れた

原文 (GenEll p.5):
> Let us refer to a compact subset of a topological space which is equal

★`IsCompactDomain`（コンパクトかつ内部の閉包に等しい）を本ファイルで定義する。
★★原文はこの語を `V^arc` 側と `V^non` 側の**両方**で使うので、
アルキメデス側からも参照できるよう一般の位相空間で書く。

## ★逸脱の記録（明示）

★★**`Gal(ℚ̄_p/ℚ_p)`-安定性と「各 `X(K)` との交わりが compact domain」は
`K_v` の側の条件として課していない。**
アルキメデス側（`CompactBound.lean` の `BoundedByArch`）が
`ι_X`-安定性やコンパクト性を `K` に課していないのと**同じ流儀**である ——
それらは利用者が課す条件であって、「囲われている」という述語の定義には要らない。
★`IsCompactDomain` は定義したので、課したいときは課せる。

★★**原文の最後の 1 文**

> the bounding domains of a compactly bounded subset ... are completely determined
> by the compactly bounded subset itself

は「well-known approximation results in elementary number theory」に依るもので、
**本ファイルは取っていない**（`.needs` に記録）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-! ## ★1. compact domain -/

/-- ★★**compact domain**（原文 p.5 の定義そのもの）。

原文 (GenEll p.5):
> Let us refer to a compact subset of a topological space which is equal

★「内部の閉包に等しいコンパクト集合」。★アルキメデス側・非アルキメデス側の両方で使う。 -/
def IsCompactDomain {T : Type*} [TopologicalSpace T] (S : Set T) : Prop :=
  IsCompact S ∧ S = closure (interior S)

theorem IsCompactDomain.isCompact {T : Type*} [TopologicalSpace T] {S : Set T}
    (h : IsCompactDomain S) : IsCompact S := h.1

theorem IsCompactDomain.eq_closure_interior {T : Type*} [TopologicalSpace T] {S : Set T}
    (h : IsCompactDomain S) : S = closure (interior S) := h.2

/-- ★compact domain は閉である（コンパクトかつ閉包に等しいから）。 -/
theorem IsCompactDomain.isClosed {T : Type*} [TopologicalSpace T] {S : Set T}
    (h : IsCompactDomain S) : IsClosed S := by
  rw [h.eq_closure_interior]; exact isClosed_closure

/-! ## ★2. `X(ℚ̄_p)` の点 -/

/-- ★**`Spec ℚ̄_p ⟶ Spec 𝓞_F`** —— 埋め込み `σ : F → ℚ̄_p` が定める射。

★アルキメデス側の `archSpecMap`（`ArchPoint.lean`）の非アルキメデス版である。

★★**索引は `F → ℚ̄_p`（`𝓞_F → ℚ̄_p` ではない）**。理由は 2 つ:
原文が「the set of **`[F:ℚ]` points**」と数えているのは体の埋め込みだからであり、
★★★体からの射は自動的に単射なので、`K/F` への**延長**が `IsAlgClosed.lift` で直に出る。 -/
noncomputable def padicSpecMap (F : Type) [Field F] [NumberField F] (p : ℕ) [Fact p.Prime]
    (σ : F →+* AlgebraicClosure ℚ_[p]) :
    Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ specRingOfIntegers F :=
  Spec.map (CommRingCat.ofHom (σ.comp (algebraMap (𝓞 F) F)))

/-- ★★**`x ∈ X(F)` と埋め込み `σ` が定める `X(ℚ̄_p)` の点**。

原文 (GenEll p.6):
> for the subset of points x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, such that for each -/
noncomputable def padicPoint {F : Type} [Field F] [NumberField F] {p : ℕ} [Fact p.Prime]
    {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X)
    (σ : F →+* AlgebraicClosure ℚ_[p]) :
    (Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ X) :=
  padicSpecMap F p σ ≫ xF

/-- ★★**`x` が定める `X(ℚ̄_p)` の点全体**。

原文の「the set of `[F:ℚ]` points of `X(ℚ̄_v)` determined by `x`」である。
★アルキメデス側の `archPointSet` と同じ形（埋め込みの像の全体）。 -/
def padicPointSet (F : Type) [Field F] [NumberField F] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    Set (Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ X) :=
  Set.range (padicPoint xF)

/-! ## ★3. 囲われていること -/

/-- ★★**`x` が `K` に囲われている**（原文「is contained in `K_v`」の非アルキメデス側）。

★★`Example 1.3, (ii)` の定義本体であり、**posit ではない**。 -/
def BoundedByNonArch (F : Type) [Field F] [NumberField F] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} (K : Set (Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ X))
    (xF : specRingOfIntegers F ⟶ X) : Prop :=
  padicPointSet F p xF ⊆ K

theorem boundedByNonArch_iff (F : Type) [Field F] [NumberField F] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} (K : Set (Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ X))
    (xF : specRingOfIntegers F ⟶ X) :
    BoundedByNonArch F p K xF
      ↔ ∀ σ : F →+* AlgebraicClosure ℚ_[p], padicPoint xF σ ∈ K := by
  constructor
  · intro h σ; exact h ⟨σ, rfl⟩
  · rintro h _ ⟨σ, rfl⟩; exact h σ

@[simp] theorem boundedByNonArch_univ (F : Type) [Field F] [NumberField F] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    BoundedByNonArch F p
      (Set.univ : Set (Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ X)) xF :=
  fun _ _ => trivial

theorem BoundedByNonArch.mono (F : Type) [Field F] [NumberField F] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} {K L : Set (Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ X)}
    (hKL : K ⊆ L) {xF : specRingOfIntegers F ⟶ X} (h : BoundedByNonArch F p K xF) :
    BoundedByNonArch F p L xF :=
  h.trans hKL

/-! ## ★4. 定義体を上げても変わらないこと（関手性の側） -/

/-- ★★埋め込みの制限に沿った関手性 —— `𝓞_K → ℚ̄_p` を `𝓞_F` へ制限する。 -/
theorem padicPoint_comp (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
    [Algebra F K] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X)
    (σ : K →+* AlgebraicClosure ℚ_[p]) :
    padicPoint (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) σ
      = padicPoint xF (σ.comp (algebraMap F K)) := by
  show padicSpecMap K p σ ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF
    = padicSpecMap F p (σ.comp (algebraMap F K)) ≫ xF
  rw [← Category.assoc]
  congr 1
  show Spec.map (CommRingCat.ofHom (σ.comp (algebraMap (𝓞 K) K)))
      ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))
    = Spec.map (CommRingCat.ofHom ((σ.comp (algebraMap F K)).comp (algebraMap (𝓞 F) F)))
  rw [← Spec.map_comp]
  -- ★タワー: `𝓞_F → 𝓞_K → K` = `𝓞_F → F → K`
  have htower : (σ.comp (algebraMap (𝓞 K) K)).comp (algebraMap (𝓞 F) (𝓞 K))
      = (σ.comp (algebraMap F K)).comp (algebraMap (𝓞 F) F) := by
    ext x
    show σ (algebraMap (𝓞 K) K (algebraMap (𝓞 F) (𝓞 K) x))
      = σ (algebraMap F K (algebraMap (𝓞 F) F x))
    rw [← IsScalarTower.algebraMap_apply, ← IsScalarTower.algebraMap_apply]
  rw [← htower]
  rfl

/-- ★★★**定義体を上げても新しい点は出ない**。

原文 (GenEll p.6):
> for the subset of points x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, such that for each

★この向きは**関手性だけ**で出る（`ℚ̄_p` が代数閉であることは要らない）。
★★逆向きは `exists_padic_extension`(代数閉性)で取ってある —— 併せて `padicPointSet_baseChange`。 -/
theorem padicPointSet_baseChange_subset (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    padicPointSet K p (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF)
      ⊆ padicPointSet F p xF := by
  rintro _ ⟨σ, rfl⟩
  exact ⟨σ.comp (algebraMap F K), (padicPoint_comp F K p xF σ).symm⟩

/-- ★★★**埋め込みは延びる** —— `ℚ̄_p` が代数閉で `K/F` が代数的だから。

★★ここが「索引を `F → ℚ̄_p` に取った」ことの効き目である ——
体からの射なので `IsAlgClosed.lift` がそのまま使える。 -/
theorem exists_padic_extension (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (p : ℕ) [Fact p.Prime]
    (σ : F →+* AlgebraicClosure ℚ_[p]) :
    ∃ τ : K →+* AlgebraicClosure ℚ_[p], τ.comp (algebraMap F K) = σ := by
  letI : Algebra F (AlgebraicClosure ℚ_[p]) := σ.toAlgebra
  haveI : Algebra.IsAlgebraic F K := inferInstance
  let ψ : K →ₐ[F] AlgebraicClosure ℚ_[p] := IsAlgClosed.lift
  exact ⟨ψ.toRingHom, by ext x; exact ψ.commutes x⟩

/-- ★★★★★**`x` が定める `X(ℚ̄_p)` の点の集合は定義体の取り方に依らない**。

原文 (GenEll p.6):
> for the subset of points x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, such that for each

★★★**これが「`K_V` が `X(ℚ̄)` の部分集合として意味を持つ」ことの中身**である。
`⊆` は関手性、`⊇` は `ℚ̄_p` の代数閉性（`exists_padic_extension`）。 -/
theorem padicPointSet_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    padicPointSet K p (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF)
      = padicPointSet F p xF := by
  refine Set.Subset.antisymm (padicPointSet_baseChange_subset F K p xF) ?_
  rintro _ ⟨σ, rfl⟩
  obtain ⟨τ, hτ⟩ := exists_padic_extension F K p σ
  exact ⟨τ, by rw [padicPoint_comp F K p xF τ, hτ]⟩

/-- ★★★★**囲われていることは定義体の取り方に依らない**（両向き）。 -/
theorem boundedByNonArch_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (p : ℕ) [Fact p.Prime]
    {X : Scheme.{0}} (C : Set (Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ X))
    (xF : specRingOfIntegers F ⟶ X) :
    BoundedByNonArch K p C (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF)
      ↔ BoundedByNonArch F p C xF := by
  unfold BoundedByNonArch
  rw [padicPointSet_baseChange F K p xF]

/-! ## ★5. 素点ごとの囲い込み領域を束ねる -/

/-- ★**非アルキメデス側の囲い込み領域** —— 素数 `p` と `X(ℚ̄_p)` の部分集合の対。

★原文の `V^non` の各 `v` に対する `K_v` である。 -/
structure PadicBound (X : Scheme.{0}) : Type 1 where
  /-- 素数 `p`（＝ `ℚ` の非アルキメデス素点）。 -/
  p : ℕ
  [isPrime : Fact p.Prime]
  /-- `K_v ⊆ X(ℚ̄_p)`。 -/
  K : Set (Spec (CommRingCat.of (AlgebraicClosure ℚ_[p])) ⟶ X)

/-- ★★`x` が `b` に囲われている。 -/
def PadicBound.Bounds {X : Scheme.{0}} (b : PadicBound X) (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X) : Prop :=
  haveI := b.isPrime
  BoundedByNonArch F b.p b.K xF

/-- ★★★★**`K_V` の条件そのもの** —— アルキメデス側と非アルキメデス側の連言。

原文 (GenEll p.6):
> for the subset of points x ∈ X(F ) ⊆ X(Q), where [F : Q] < ∞, such that for each

★★原文の `V` は `V(ℚ)^arc` を必ず含む有限集合なので、
`ℚ` の唯一のアルキメデス素点に対する `Karc` と、有限個の `PadicBound` の族で書ける。 -/
def BoundedByV {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (Karc : Set (complexPoints X)) (Vnon : List (PadicBound X))
    (xF : specRingOfIntegers F ⟶ X) : Prop :=
  BoundedByArch F Karc xF ∧ ∀ b ∈ Vnon, b.Bounds F xF

theorem boundedByV_of_arch {X : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (Karc : Set (complexPoints X)) (xF : specRingOfIntegers F ⟶ X)
    (h : BoundedByArch F Karc xF) : BoundedByV F Karc [] xF :=
  ⟨h, by simp⟩

/-! ### ★出典の紐付け(`.src`) -/

def IsCompactDomain.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Example 1.3, (ii)(compact domain の定義)",
    sectionId := "genell-ex-1-3" }

def padicPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Example 1.3, (ii)(x が定める X(ℚ̄_p) の点)",
    sectionId := "genell-ex-1-3" }

def BoundedByNonArch.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Example 1.3, (ii)(非アルキメデス側の囲い込み)",
    sectionId := "genell-ex-1-3" }

def BoundedByV.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Example 1.3, (ii)(K_V の条件——アルキメデス側と非アルキメデス側の連言)",
    sectionId := "genell-ex-1-3" }

def BoundedByV.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "BoundedByArch(アルキメデス側。CompactBound.lean)"
      (.inProject "ABC3" "ABC3.Found.GenEll.BoundedByArch") 6,
    .implicitStep
      ("★逸脱: K_v の側の条件(Gal(ℚ̄_p/ℚ_p)-安定性、各 X(K) との交わりが compact domain、" ++
       "ι_X-安定性、K_v ≠ 全体)は課していない。アルキメデス側と同じ流儀で、" ++
       "それらは利用者が課す条件であり「囲われている」という述語の定義には要らない。" ++
       "★IsCompactDomain は定義したので課したいときは課せる") 6,
    .folklore
      ("原文「by applying well-known approximation results in elementary number theory, " ++
       "it follows immediately that the bounding domains ... are completely determined by " ++
       "the compactly bounded subset itself」——★本ファイルは取っていない") 6 ]

end ABC3.Found.GenEll
