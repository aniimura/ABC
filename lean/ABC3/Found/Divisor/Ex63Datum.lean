/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithSurj
import ABC3.Found.Divisor.ArithPhiPerf
import ABC3.Found.Divisor.ArithFrobenioid
import ABC3.Found.FrdI.Sec6GaloisCat
import ABC3.Found.FrdI.MonoidTransport

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

variable {F Kbar}

/-! ## ★★★★★5. `Φ` は non-dilating —— `Theorem 6.4, (i)` の一条(2026-08-21)

原文は「数体の自己同型で全素点を固定するものは恒等」を `clearly` で畳む。
★★その non-dilating の側の中身は次の 3 つに割れた:

| 段 | 宣言 |
|---|---|
| `𝒟` の自己射は同型 | `finSubEndo_bijective`(`AlgHom.bijective`、有限次元) |
| 同型に沿った引き戻しは係数を変えない | `arithExtend_apply'`(`localDeg = 1`) |
| primary 元で生成され、primary 元は台が 1 点 | `closure_primary_effR_eq_top` / `effR_single_supported_of_primary` |

★★★**アルキメデス素点があるので抽象論では出ない** ——
成分が `ℝ≥0` そのものなので `c ↦ 2c` という加法的全単射が存在する。
「係数 1」は **`localDeg = 1`**、すなわち数論の入力である。 -/

omit [NumberField F] in
/-- ★**`FinSub` の自己射は全単射**(有限次元だから)。 -/
theorem finSubEndo_bijective {L : FinSub F Kbar} (σ : L ⟶ L) :
    Function.Bijective (FinSub.hom σ) := by
  haveI := L.fin
  exact AlgHom.bijective (FinSub.hom σ)

omit [NumberField F] in
theorem algOfHom_bijective {L : FinSub F Kbar} (σ : L ⟶ L) :
    Function.Bijective (@algebraMap L.toIF L.toIF _ _ (algOfHom σ)) :=
  finSubEndo_bijective σ

/-- ★★**自己射に沿った引き戻しは「素点を引き戻す」だけ**。 -/
theorem pullOf_apply {L : FinSub F Kbar} (σ : L ⟶ L) (d : ArithPlace L.toIF →₀ ℝ)
    (V : ArithPlace L.toIF) :
    pullOf σ d V = d (@resPlace L.toIF L.toIF _ _ _ _ (algOfHom σ) V) := by
  letI := algOfHom σ
  exact arithExtend_apply' (algOfHom_bijective σ) d V

/-- ★★★★★★**[FrdI] Theorem 6.4, (i)** —— **`Φ` は non-dilating**。

原文 (FrdI p.114):
> group-like monoid on Di given by the multiplicative group of the field extension

★★primary な元は台が 1 点(`effR_single_supported_of_primary`)なので、
non-dilating の仮定「`α^char(a) ≼ a`」は
**`a` の台の素点が `α` で動かない**ことを言う。
★そこで `α` が単射であることから `resPlace v = v` が出て、
`localDeg = 1` により係数も変わらない。 -/
theorem arithDatumGalois_isNonDilatingOn :
    MonoidOn.IsNonDilatingOn ((arithDatumGalois F Kbar).phi finSubOp_isOfFSMType) := by
  intro A e
  set Δ := arithDatumGalois F Kbar with hΔ
  have hmap : (Δ.phi finSubOp_isOfFSMType).map e = Δ.mapHom e := rfl
  rw [hmap]
  refine isNonDilating_of_primary_sharp (isSharp_effR _) _
    (closure_primary_effR_eq_top (Δ.coord A)) ?_
  intro a ha hprec
  obtain ⟨v, hv⟩ := effR_single_supported_of_primary (Δ.coord A) ha
  set R : Δ.primes A → Δ.primes A :=
    fun V => @resPlace _ _ _ _ _ _ (algOfHom e.unop) V with hR
  have happ : ∀ (x : effR (Δ.grp A)) (V : Δ.primes A),
      ((Δ.mapHom e x : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) V
        = ((x : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) (R V) :=
    fun x V => pullOf_apply e.unop _ V
  obtain ⟨n, hn, c, hc⟩ := hprec
  have hcoe : ((Δ.mapHom e a : effR (Δ.grp A)) : Δ.primes A →₀ ℝ)
      + ((c : effR (Δ.grp A)) : Δ.primes A →₀ ℝ)
      = n • ((a : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) :=
    congrArg Subtype.val hc
  have hzero : ∀ V : Δ.primes A, V ≠ v →
      ((a : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) (R V) = 0 := by
    intro V hV
    have h1 := congrArg (fun f : Δ.primes A →₀ ℝ => f V) hcoe
    simp only [Finsupp.add_apply, Finsupp.smul_apply] at h1
    rw [hv, Finsupp.single_eq_of_ne hV, smul_zero] at h1
    rw [← happ a V]
    have h2 := (mem_effR.mp (Δ.mapHom e a).2).2 V
    have h3 := (mem_effR.mp c.2).2 V
    linarith
  have hane : (Δ.mapHom e a) ≠ 0 := fun h =>
    ha.1 (Δ.mapHom_injective e (by rw [h, map_zero]))
  have hex : ∃ V : Δ.primes A, ((Δ.mapHom e a : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) V ≠ 0 := by
    by_contra hcon
    push_neg at hcon
    exact hane (Subtype.ext (Finsupp.ext hcon))
  obtain ⟨V₀, hV₀⟩ := hex
  have hV₀v : V₀ = v := by
    by_contra hne
    exact hV₀ ((happ a V₀).trans (hzero V₀ hne))
  rw [hV₀v] at hV₀
  have hRv : R v = v := by
    by_contra hne
    refine hV₀ ((happ a v).trans ?_)
    rw [hv]
    exact Finsupp.single_eq_of_ne hne
  refine Subtype.ext (Finsupp.ext fun V => ?_)
  rcases eq_or_ne V v with rfl | hV
  · rw [happ a V, hRv]
  · rw [happ a V, hzero V hV, hv]
    exact (Finsupp.single_eq_of_ne hV).symm

/-- ★★★★★locator —— `Theorem 6.4, (i)` の non-dilating。 -/
def arithDatumGalois_isNonDilatingOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Φ は non-dilating(自己同型は素点を係数 1 で移す)",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.Divisor
