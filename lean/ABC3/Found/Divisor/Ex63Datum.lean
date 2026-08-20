/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithSurj
import ABC3.Found.Divisor.ArithPhiPerf
import ABC3.Found.Divisor.ArithFrobenioid
import ABC3.Found.FrdI.Sec6GaloisCat

/-!
# `Example 6.3` の `ArithDatum` を `𝒟 = B(G)⁰` の上で実現する

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★底の圏は `Sec6GaloisCat.lean` の `(FinSub F F̄)ᵒᵖ`

原文の `𝒟 = B(G)⁰` を、`Sec6GaloisCat.lean` は
「`F̄/F` の有限部分拡大を対象、`F`-代数射を射とする圏の**反対圏**」として実現した。
★その圏について **連結**(`finSubOp_isConnected`)・**totally epimorphic**
(`finSubOp_totallyEpimorphic`)・**of FSM-type**(`finSubOp_isOfFSMType`)は
すでに在庫にある。

## ★★★底を `ℚ` に取ってはいけない(実測 2026-08-21)

★★`F = ℚ` に取ると **`Algebra ℚ (IntermediateField ℚ F̄)` に diamond が立つ** ——
`DivisionRing.toRatAlgebra` と `IntermediateField.algebra'` が
**定義的に等しくない**ので、`AlgHom ℚ` の項が型検査を通らない。
★**底を一般の数体 `F` に取れば diamond は消える**(原文も `F` は一般の数体である)。

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `nfFinSub` | `L ∈ Ob(B(G)⁰)` は数体(`[L:ℚ] = [L:F]·[F:ℚ]`) |
| `algOfHom` | 射 `f : L → M` が定める `Algebra L M` |
| `pullOf` | `f` に沿った算術因子の引き戻し |
| `pullOf_id` / `pullOf_comp` | `ArithTower.lean` の `arithExtend_id` / `arithExtend_comp` |
| `arithDatumGalois` | ★★★★**`ArithDatum (FinSub F F̄)ᵒᵖ`** |

★★残るのは**有理関数の単系 `B(L) = L^×`** の側(`bmon` と `divB`)だけである。
-/

namespace ABC3.Found.Divisor

open CategoryTheory NumberField ABC3.Found.FrdI

open scoped NumberField

variable {F Kbar : Type} [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]

/-! ## ★1. 対象は数体 -/

/-- ★**`B(G)⁰` の対象は数体** —— `[L:ℚ] = [L:F]·[F:ℚ] < ∞`。 -/
noncomputable instance nfFinSub (L : FinSub F Kbar) : NumberField L.toIF := by
  haveI := L.fin
  haveI : FiniteDimensional ℚ L.toIF := Module.Finite.trans F L.toIF
  exact ⟨⟩

/-! ## ★2. 射が定める代数構造 -/

/-- ★**射 `f : L → M` が定める `Algebra L M`**。

★★`f = 𝟙` のとき `Algebra.id` に**定義的に等しい**のが要点である
(`RingHom.toAlgebra (RingHom.id _)` がまさに `Algebra.id`)。 -/
@[reducible] noncomputable def algOfHom {L M : FinSub F Kbar} (f : L ⟶ M) :
    Algebra L.toIF M.toIF :=
  (FinSub.hom f).toRingHom.toAlgebra

omit [NumberField F] in
theorem algOfHom_id (L : FinSub F Kbar) : algOfHom (𝟙 L) = Algebra.id L.toIF := rfl

/-! ## ★3. 引き戻し -/

/-- ★★**算術因子の引き戻し** —— `ArithFunctor.lean` の `arithExtend`。 -/
noncomputable def pullOf {L M : FinSub F Kbar} (f : L ⟶ M) :
    (ArithPlace L.toIF →₀ ℝ) →+ (ArithPlace M.toIF →₀ ℝ) :=
  letI := algOfHom f
  arithExtend

/-- ★**恒等射での引き戻しは恒等**(`ArithTower.lean` の `arithExtend_id`)。 -/
theorem pullOf_id (L : FinSub F Kbar) (x : ArithPlace L.toIF →₀ ℝ) :
    pullOf (𝟙 L) x = x :=
  arithExtend_id x

/-- ★★**引き戻しは合成と両立する**(`ArithTower.lean` の `arithExtend_comp`)。 -/
theorem pullOf_comp {L M N : FinSub F Kbar} (f : L ⟶ M) (g : M ⟶ N)
    (x : ArithPlace L.toIF →₀ ℝ) :
    pullOf (f ≫ g) x = pullOf g (pullOf f x) := by
  letI := algOfHom f
  letI := algOfHom g
  letI := algOfHom (f ≫ g)
  haveI : IsScalarTower L.toIF M.toIF N.toIF :=
    IsScalarTower.of_algebraMap_eq (fun _ => rfl)
  exact arithExtend_comp x

/-! ## ★4. `ArithDatum` -/

variable (F Kbar)

/-- ★★★★★★**`Example 6.3` の算術因子のデータ**を `𝒟 = B(G)⁰` の上で実現したもの。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model] -/
noncomputable def arithDatumGalois : ArithDatum.{0, 0, 0} (FinSub F Kbar)ᵒᵖ where
  primes A := ArithPlace A.unop.toIF
  pull {_ _} α := pullOf α.unop
  pull_id A x := pullOf_id A.unop x
  pull_comp {_ _ _} α β x := pullOf_comp α.unop β.unop x
  grp A := arithDivGroup A.unop.toIF
  pull_mem {_ _} α {_} hx := by
    letI := algOfHom α.unop
    exact arithExtend_mem_arithDivGroup hx
  pull_nonneg {_ _} α {_} hx := by
    letI := algOfHom α.unop
    exact arithExtend_nonneg hx
  pull_inj {_ _} α := by
    letI := algOfHom α.unop
    exact arithExtend_injective
  gen A := isGenSubgroupR_arithDivGroup
  locMono A := isLocallyMonoprimeR_arithDivGroup
  coord A := isCoordwiseR_arithDivGroup

variable [IsGalois F Kbar]

/-- ★★★**因子の単系 `Φ` が `𝒟` 上の monoid になる**(`Theorem 5.2` の入力の半分)。 -/
noncomputable def phiGalois : MonoidOn.{0, 0, 0} (FinSub F Kbar)ᵒᵖ :=
  (arithDatumGalois F Kbar).phi finSubOp_isOfFSMType

/-- ★★**`Φ(L)` は divisorial**。 -/
theorem phiGalois_isDivisorialOn : (phiGalois F Kbar).IsDivisorialOn :=
  (arithDatumGalois F Kbar).phi_isDivisorialOn finSubOp_isOfFSMType

/-- ★★★**`Φ(L)` は perf-factorial**(`Example 6.3` の `one verifies immediately`)。 -/
theorem phiGalois_isPerfFactorial (A : (FinSub F Kbar)ᵒᵖ) :
    IsPerfFactorial ((phiGalois F Kbar).val A) :=
  (arithDatumGalois F Kbar).phi_isPerfFactorial finSubOp_isOfFSMType A

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Example 6.3` の算術因子のデータを `B(G)⁰` の上で実現する段。 -/
def arithDatumGalois.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 算術因子のデータを 𝒟 = B(G)⁰ の上で実現する",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
