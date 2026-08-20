/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.Divisor.SchemeWeil
import ABC3.Found.Divisor.CartierMonoid

/-!
# 因子論の第 2 層 —— Cartier 因子と `Q`-Cartier 性(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109(`Example 6.1`)。

原文 (FrdI p.109):
> DK a set of Q-Cartier prime divisors on V . The connected objects of the Galois category

原文 (FrdI p.109):
> for the monoid of Cartier effective divisors D on V [L] with support in DL [i.e., D

## ★★なぜ原文が Cartier に限るのか

★**Weil 因子は一般には引き戻せない**が、**Cartier 因子は引き戻せる**。
`Example 6.1` は `L ↦ Φ(L)` を `𝒟 = B(G)⁰` 上の**単系(関手)**にする必要があるので、
引き戻しが要る —— そこが Cartier に限る理由である
(`Theorem 6.2, (i)` の仮定 (a) も同じ所で効く)。

## ★★単系論の側はすでに全部閉じている

`Found/Divisor/Cartier*.lean` の 5 ファイルで、`S` を素因子の型、
`Γ ≤ ℤ[S]` を Cartier 因子の群、`Φ := Γ ∩ ℤ≥0[S]` と置いたときの

* `isDivisorial_effSub` —— ★**`Q`-Cartier 性は要らない**(効くのは `Γ` が部分群だけ)
* `effSubGpEquiv : Gp (effSub Γ) ≃+ Γ`
* `mprec_effSub_iff` —— `a ⪯ b ⇔ supp a ⊆ supp b`
* `effSubPrimeEquiv : Prime (effSub Γ) ≃ S`、`pfEquivNonneg : Pf (effSub Γ) ≃+ ℚ≥0[S]`
* `isPerfFactorial_effSub` —— perf-factorial(11 条全部)

を証明してある。★**したがってこの層に残るのは「幾何の `Γ` が実際に部分群であり
`Q`-Cartier である」ことだけ**である(鎖 `cartier` の `cartier-subgroup` / `qcartier`)。

## ★在庫

`Found/Arakelov/PicDivisor.lean` に `IsCartier X D`(`X.IdealSheafData` 版)が在る。
★本ファイルの `IsCartierDiv`(Weil 因子版)との橋は未実装である。
-/

namespace ABC3.Skeleton.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Meta

universe u

variable (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]

/-! ## ★1. Cartier 因子(鎖 `cartier` の `cartier-def`) -/

/-- ★**Cartier 因子** —— 局所的に有理函数の因子であるもの。

★各点のまわりで、ある `f ∈ K(X)^×` があって `D` がその近傍で `div(f)` に一致する。 -/
def IsCartierDiv (D : WeilDiv X) : Prop :=
  ∀ x : X, ∃ (U : X.Opens) (_ : x ∈ U) (f : (X.functionField)ˣ),
    ∀ y : PrimeDivisorPt X, y.1 ∈ U → D y = ordAtDiv X y (f : X.functionField)

/-- ★**`Q`-Cartier 因子** —— ある正の倍数が Cartier。 -/
def IsQCartierDiv (D : WeilDiv X) : Prop :=
  ∃ n : ℕ, 0 < n ∧ IsCartierDiv X (n • D)

/-! ## ★2. Cartier 因子は部分群をなす(鎖 `cartier` の `cartier-subgroup`)

★★**ここが単系論への入口である。**`Found/Divisor/CartierMonoid.lean` は
`Γ : AddSubgroup (S →₀ ℤ)` しか要求しないので、この 3 条が出れば
`Example 6.1` の単系論(divisorial・`Φ^gp ≃ Γ`・`Prime ≃ D_L`・perf-factorial)は
**すべて自動で付く**。 -/

/-- ★`0` は Cartier。 -/
theorem isCartierDiv_zero : IsCartierDiv X 0 := by
  sorry

/-- ★Cartier 因子の和は Cartier。 -/
theorem isCartierDiv_add {D E : WeilDiv X} (_hD : IsCartierDiv X D) (_hE : IsCartierDiv X E) :
    IsCartierDiv X (D + E) := by
  sorry

/-- ★Cartier 因子の逆元は Cartier。 -/
theorem isCartierDiv_neg {D : WeilDiv X} (_hD : IsCartierDiv X D) :
    IsCartierDiv X (-D) := by
  sorry

/-- ★★**Cartier 因子の群** —— `Example 6.1` の `Φ(L)^gp ⊆ ℤ[D_L]` の中身。 -/
def cartierSubgroup : AddSubgroup (WeilDiv X) where
  carrier := {D | IsCartierDiv X D}
  zero_mem' := isCartierDiv_zero X
  add_mem' hD hE := isCartierDiv_add X hD hE
  neg_mem' hD := isCartierDiv_neg X hD

/-! ## ★3. `Q`-Cartier 性を単系論の語彙へ(鎖 `cartier` の `qcartier`)

★`Found/Divisor/CartierMonoid.lean` の `IsQCartierSubgroup Γ := ∀ s, ∃ n>0, single s n ∈ Γ`
と、幾何側の「各素因子が `Q`-Cartier」が一致する、というのがこの節点である。 -/

/-- ★★**各素因子が `Q`-Cartier ⟹ `Γ` は `Q`-Cartier 部分群**。

★`Γ` は `S`(素因子の部分集合)へ制限した Cartier 因子の群。
★これが `Found/Divisor/Cartier*.lean` の全定理の前提 `hQ` を与える。 -/
theorem isQCartierSubgroup_of_forall_isQCartier
    {S : Type u} (ι : S → PrimeDivisorPt X) (_hι : Function.Injective ι)
    (Γ : AddSubgroup (S →₀ ℤ))
    (_hΓ : ∀ x : S →₀ ℤ, x ∈ Γ ↔ IsCartierDiv X (Finsupp.mapDomain ι x))
    (_hQ : ∀ s : S, IsQCartierDiv X (Finsupp.single (ι s) 1)) :
    ABC3.Found.FrdI.IsQCartierSubgroup Γ := by
  sorry

/-! ## ★4. Cartier 因子の引き戻し(鎖 `cartier` の `cartier-pullback`)

★★**`Example 6.1` の関手性の本体**である。 -/

/-- ★★**支配射に沿った Cartier 因子の引き戻し**。

★`Φ(L) → Φ(M)`(`L ⊆ M`)と `Theorem 6.2, (i)` の `Φ₁ → Φ₂|𝒟₁` は、
どちらもこれである。 -/
noncomputable def pullbackCartier {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (_ψ : Y ⟶ X)
    (_D : WeilDiv X) (_hD : IsCartierDiv X _D) : WeilDiv Y :=
  sorry

/-- ★引き戻しは加法的。 -/
theorem pullbackCartier_add {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X)
    {D E : WeilDiv X} (hD : IsCartierDiv X D) (hE : IsCartierDiv X E) :
    pullbackCartier X ψ (D + E) (isCartierDiv_add X hD hE)
      = pullbackCartier X ψ D hD + pullbackCartier X ψ E hE := by
  sorry

/-- ★引き戻しは Cartier 性を保つ。 -/
theorem isCartierDiv_pullbackCartier {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X)
    {D : WeilDiv X} (hD : IsCartierDiv X D) :
    IsCartierDiv Y (pullbackCartier X ψ D hD) := by
  sorry

/-- ★引き戻しは有効性を保つ(`Φ` を `Φ` へ移すために要る)。 -/
theorem pullbackCartier_nonneg {Y : Scheme.{u}} [IsIntegral Y]
    [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X)
    {D : WeilDiv X} (hD : IsCartierDiv X D) (_hpos : 0 ≤ D) :
    0 ≤ pullbackCartier X ψ D hD := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def IsCartierDiv.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Cartier divisor",
    sectionId := "frdi-example-6-1" }

def IsQCartierDiv.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Q-Cartier prime divisor",
    sectionId := "frdi-example-6-1" }

def isCartierDiv_zero.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Cartier 因子の群(零元)",
    sectionId := "frdi-example-6-1" }

def isCartierDiv_zero.needs : List ProofObligation :=
  [ .derivation "`f = 1` を取れば `div(1) = 0`" 109 ]

def isCartierDiv_add.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Cartier 因子の群(和)",
    sectionId := "frdi-example-6-1" }

def isCartierDiv_add.needs : List ProofObligation :=
  [ .derivation "2 つの近傍の共通部分を取り、`f·g` を当てる(`ordAtDiv_mul`)" 109,
    .citation "[ABC3]" "ordAtDiv_mul" (.inProject "ABC3" "ABC3.Skeleton.Divisor.ordAtDiv_mul") 109 ]

def isCartierDiv_neg.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Cartier 因子の群(逆元)",
    sectionId := "frdi-example-6-1" }

def isCartierDiv_neg.needs : List ProofObligation :=
  [ .derivation "`f⁻¹` を当てる" 109 ]

def cartierSubgroup.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Φ(L)^gp ⊆ ℤ[D_L] は Cartier 因子の群",
    sectionId := "frdi-example-6-1" }

def isQCartierSubgroup_of_forall_isQCartier.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — K-Q-Cartier ⟹ IsQCartierSubgroup",
    sectionId := "frdi-example-6-1" }

/-- ★★これが出れば `Found/Divisor/Cartier*.lean` の**全定理の前提が揃う**。 -/
def isQCartierSubgroup_of_forall_isQCartier.needs : List ProofObligation :=
  [ .citation "[ABC3]" "IsQCartierSubgroup(単系側の定義)"
      (.inProject "ABC3" "ABC3.Found.FrdI.IsQCartierSubgroup") 109,
    .derivation "`n·[s]` が Cartier なら `Finsupp.single s n ∈ Γ`(`mapDomain` の単射性で戻す)" 109,
    .implicitStep "★原文は `K-Q-Cartier` を**仮定**として置く(「we shall assume that DK is K-Q-Cartier」)" 109 ]

def pullbackCartier.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — pulling back divisors",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_add.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しの加法性",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_add.needs : List ProofObligation :=
  [ .derivation "局所的に `f ↦ f ∘ ψ` を取るだけなので、積が和に移る" 110 ]

def isCartierDiv_pullbackCartier.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しは Cartier 性を保つ",
    sectionId := "frdi-thm-6-2" }

def isCartierDiv_pullbackCartier.needs : List ProofObligation :=
  [ .derivation "局所主因子の引き戻しは局所主因子" 110,
    .implicitStep
      "★原文は「by pulling back [Cartier] divisors and rational functions via ψ, we obtain compatible natural transformations」で畳む" 110 ]

def pullbackCartier_nonneg.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2, (i) — 引き戻しは有効性を保つ",
    sectionId := "frdi-thm-6-2" }

def pullbackCartier_nonneg.needs : List ProofObligation :=
  [ .derivation "支配射なら局所方程式が 0 にならないので、位数は非負のまま" 110 ]

end ABC3.Skeleton.Divisor
