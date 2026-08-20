import ABC3.Found.GaloisRep.TranslateComp

/-!
# Galois (G5) 第 122 ブロック —— **★★★★★★2 等分点でも単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★分解を作って合成則に流す

第 121 で「分解 `Q₃ = Q₁ + Q₂` があれば `τ_{Q₃}` は単射」を得た。
★その**分解を実際に作る**:

    Q₁ := 2 等分点でないアフィン点(1 つあればよい)
    Q₂ := Q₃ − Q₁

★★`2Q₂ = 2Q₃ − 2Q₁ = −2Q₁ ≠ 0` なので **`Q₂` は自動的に 2 等分点でない**。
★★★`Q₂ ≠ 0` もそこから出る。

## ★★★★これで平行移動の単射性が完成した

    2 等分点でない Q  → 第 117(`−Q` での 1 点評価)
    2 等分点の Q     → 第 121 + 本ブロック(分解して合成)

★残るのは「2 等分点でないアフィン点が 1 つ存在すること」だけであり、
`E[2]` は高々 4 点(第 65-72 の `E[n] ≃ (ℤ/n)²`)なので体が十分大きければ取れる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `two_nsmul_eq_zero_iff_negY` | ★★2 等分点 ⟺ `negY x y = y` |
| `exists_decomp_of_twoTorsion` | ★★★★★**2 等分点の分解** |
| `toFF` / `toFF_some` | ★点を関数体上へ送る群準同型 |
| `translateHom_injective_all` | ★★★★★★**つねに単射** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

/-! ## ★★2 等分点の判定(`negY` の形) -/

/-- ★★**点が 2 等分点であること ⟺ `negY x y = y`**。

★第 31 ブロックの `two_smul_eq_zero_iff`(`Ψ₂Sq` の形)と同値だが、
本ブロックでは分解の議論に `negY` の形が要る。 -/
theorem two_nsmul_eq_zero_iff_negY {K : Type} [Field K] [DecidableEq K]
    {V : WeierstrassCurve.Affine K} {x y : K} (h : V.Nonsingular x y) :
    (2 : ℕ) • (Point.some x y h) = 0 ↔ V.negY x y = y := by
  rw [two_smul, add_eq_zero_iff_eq_neg, WeierstrassCurve.Affine.Point.neg_some,
    WeierstrassCurve.Affine.Point.some.injEq]
  constructor
  · rintro ⟨-, hy⟩
    exact hy.symm
  · intro hy
    exact ⟨rfl, hy.symm⟩

/-! ## ★★★★★2 等分点の分解 -/

/-- ★★★★★**2 等分点は 2 等分点でない 2 点の和に分解できる**
(2 等分点でないアフィン点が 1 つあれば)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Q₂ := Q₃ − Q₁` と置けば `2Q₂ = −2Q₁ ≠ 0` なので自動的に 2 等分点でない。 -/
theorem exists_decomp_of_twoTorsion {K : Type} [Field K] [DecidableEq K]
    {V : WeierstrassCurve.Affine K} {x₃ y₃ : K} (h₃ : V.Nonsingular x₃ y₃)
    (hT : V.negY x₃ y₃ = y₃)
    (hex : ∃ (x₁ y₁ : K) (_ : V.Nonsingular x₁ y₁), V.negY x₁ y₁ ≠ y₁) :
    ∃ (x₁ y₁ x₂ y₂ : K) (h₁ : V.Nonsingular x₁ y₁) (h₂ : V.Nonsingular x₂ y₂),
      V.negY x₁ y₁ ≠ y₁ ∧ V.negY x₂ y₂ ≠ y₂ ∧
      Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂ = Point.some x₃ y₃ h₃ := by
  obtain ⟨x₁, y₁, h₁, k₁⟩ := hex
  set P₁ := Point.some x₁ y₁ h₁ with hP₁
  set P₃ := Point.some x₃ y₃ h₃ with hP₃
  have h2P₁ : (2 : ℕ) • P₁ ≠ 0 := fun hc => k₁ ((two_nsmul_eq_zero_iff_negY h₁).1 hc)
  have h2P₃ : (2 : ℕ) • P₃ = 0 := (two_nsmul_eq_zero_iff_negY h₃).2 hT
  have h2P₂ : (2 : ℕ) • (P₃ - P₁) ≠ 0 := by
    intro hc
    refine h2P₁ ?_
    rw [smul_sub, h2P₃, zero_sub, neg_eq_zero] at hc
    exact hc
  have hP₂ne : P₃ - P₁ ≠ 0 := by
    intro hc
    refine h2P₂ ?_
    rw [hc, smul_zero]
  obtain ⟨x₂, y₂, h₂, hP₂eq⟩ := exists_coords_of_ne_zero (P₃ - P₁) hP₂ne
  refine ⟨x₁, y₁, x₂, y₂, h₁, h₂, k₁, ?_, ?_⟩
  · intro hc
    refine h2P₂ ?_
    rw [hP₂eq]
    exact (two_nsmul_eq_zero_iff_negY h₂).2 hc
  · rw [← hP₂eq]
    abel

/-! ## ★★★★★★つねに単射 -/

variable {F : Type} [Field F] [DecidableEq F]

/-- ★点を関数体上へ送る群準同型。 -/
noncomputable def toFF (W : WeierstrassCurve.Affine F) :
    W.Point →+ (W.map (algebraMap F W.FunctionField)).Point :=
  WeierstrassCurve.Affine.Point.map (W' := W) (Algebra.ofId F W.FunctionField)

theorem toFF_some (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) :
    toFF W (Point.some x y h)
      = Point.some (algebraMap F W.FunctionField x) (algebraMap F W.FunctionField y)
        (mapNonsingular W h) :=
  WeierstrassCurve.Affine.Point.map_some _ _

/-- ★★★★★★**平行移動の環準同型はつねに単射である**
(2 等分点でないアフィン点が 1 つあれば)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★2 等分点でない場合は第 117(`−Q` での 1 点評価)、
★★2 等分点の場合は分解して第 121 の合成則に流す。 -/
theorem translateHom_injective_all (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀)
    (hex : ∃ (x₁ y₁ : F) (_ : W.Nonsingular x₁ y₁), W.negY x₁ y₁ ≠ y₁) :
    Function.Injective (translateHom W hQ) := by
  by_cases h2 : W.negY x₀ y₀ = y₀
  · obtain ⟨x₁, y₁, x₂, y₂, h₁, h₂, k₁, k₂, hsum⟩ := exists_decomp_of_twoTorsion hQ h2 hex
    refine translateHom_injective_of_decomp W h₁ h₂ hQ k₁ k₂ ?_
    have hmap := congrArg (toFF W) hsum
    rw [map_add, toFF_some, toFF_some, toFF_some] at hmap
    exact hmap
  · exact translateHom_injective W hQ h2

/-! ## ★出典の紐付け(`.src`) -/

def translateHom_injective_all.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の単射性が 2 等分点でも成り立つこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
