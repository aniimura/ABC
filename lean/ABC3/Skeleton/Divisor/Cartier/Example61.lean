/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.Divisor.SchemeWeil
import ABC3.Found.Divisor.CartierMonoid

/-!
# Cartier —— `[FrdI] Example 6.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
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

end ABC3.Skeleton.Divisor
