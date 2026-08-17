import ABC3.Found.GenEll.CartierPullback
import ABC3.Found.GenEll.DegMul

/-!
# [GenEll] Definition 1.1 の足場 —— **ℂ-点**と Green 関数の引き戻し(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

## ★★★算術直線束の引き戻しは 2 つの成分に分かれる

`x_F^* L̄` は `Spec 𝓞_F` 上の算術直線束、すなわち `ADiv F` の元である。
`ADiv F = (FinitePlace F →₀ ℤ) × (InfinitePlace F →₀ ℝ)` なので、
引き戻しも **2 つの成分**に分かれる:

| 成分 | 由来 | 状態 |
|---|---|---|
| 有限素点側 | Cartier 因子 `D` の引き戻し | ★`CartierPullback.lean` の `pullbackIdeal` |
| 無限素点側 | **Green 関数を `x` の ℂ-点で評価** | ★★**本ファイル** |

★★**無限素点側は完全に構成できる。** 解析空間も GAGA も要らない——
要るのは「`x ∈ X(F)` と無限素点 `v` から ℂ-点を作る」ことだけであり、
それは **`𝓞_F → ℂ` の `Spec` を取って合成する**だけである。

## ★★★加法性は無限素点側では**無条件に**成り立つ

Green 関数はテンソル積で**足される**(`|·|_{L⊗M} = |·|_L · |·|_M` の対数)。
評価は線形なので、`archADiv` は `g` について加法的である。
★★これは仮定を 1 つも要しない——`archADiv_add` を参照。

★有限素点側は `IdealSheafData.comap` の積の保存が要り、そこはまだ穴である
(`DegMul.lean` の冒頭に記録)。

## ★mathlib 実測(2026-08-17)

`NumberField.InfinitePlace.embedding : InfinitePlace K → (K →+* ℂ)` は在る。
`InfinitePlace K` は `Fintype` である。
★**スキームの ℂ-点(`Spec ℂ ⟶ X`)を数論的な埋め込みから作る補題は無い**——本ファイルで作る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★ℂ-点 -/

/-- ★**スキームの ℂ-点** `X(ℂ)`。

★原文の `X^arc` の**台集合**である。位相も正則構造も入れていない——
`ht` の定義に要らないからである(`ProjTopology.lean` の測定を参照)。 -/
abbrev complexPoints (X : Scheme.{0}) : Type _ := Spec (CommRingCat.of ℂ) ⟶ X

/-- ★無限素点 `v` が定める `𝓞_F → ℂ`。

`v.embedding : F →+* ℂ` に `𝓞_F → F` を前合成する。 -/
noncomputable def archRingHom (v : InfinitePlace F) : (𝓞 F) →+* ℂ :=
  v.embedding.comp (algebraMap (𝓞 F) F)

/-- ★**`Spec ℂ ⟶ Spec 𝓞_F`** —— 無限素点が定める射。 -/
noncomputable def archSpecMap (v : InfinitePlace F) :
    Spec (CommRingCat.of ℂ) ⟶ specRingOfIntegers F :=
  Spec.map (CommRingCat.ofHom (archRingHom F v))

variable {F}

/-- ★★**`x ∈ X(F)` と無限素点 `v` が定める ℂ-点**。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★原文の「the set of [F : Q] points of X^arc determined by x」の 1 つ 1 つである
(物理 p.6、`Example 1.3, (ii)` でも同じものが使われる)。 -/
noncomputable def archPoint {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X)
    (v : InfinitePlace F) : complexPoints X :=
  archSpecMap F v ≫ xF

/-! ## ★★Green 関数とその引き戻し -/

/-- ★**Green 関数** —— `X(ℂ)` 上の実数値関数。

★原文の hermitian 計量 `|−|_L` の `−log` にあたる。
★★**計量そのものではなく対数を取った形で持つ**——
テンソル積が積でなく**和**になり、`ADiv` の加法群構造とそのまま合うからである。 -/
abbrev GreenFn (X : Scheme.{0}) : Type _ := complexPoints X → ℝ

variable (F)

open scoped Classical in
/-- ★★**Green 関数の引き戻し** —— 各無限素点で `x` の ℂ-点における値を取り、
**局所次数 `mult v` で重みを付ける**。

★`InfinitePlace F` は有限なので `Finsupp` に載る。

## ★★★`mult` の重みが要る理由(2026-08-17 深夜に訂正)

★当初 `g (archPoint xF v)` だけを係数にしていた。**誤りだった。**

`ADiv` の底変換(`BaseChange.lean` の `baseChangeArc`)は
係数に `mult w / mult v` を掛ける。★★`archPoint` の値は底変換で**変わらない**ので、
`mult` の重みが無いと `archADiv` は `baseChange` の形にならず、
**高さが定義体の取り方に依ってしまう**。

★★★`Σ_v mult v = [F : ℚ]` なので、`mult` を重みに入れて `[F:ℚ]` で割ると
**重み付き平均**になり、底変換で不変になる。

★**底変換を証明しようとして初めて露見した誤りである。** -/
noncomputable def archADiv {X : Scheme.{0}} (g : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) : InfinitePlace F →₀ ℝ :=
  Finsupp.onFinset Finset.univ
    (fun v => (InfinitePlace.mult v : ℝ) * g (archPoint xF v)) (by
    intro v _; exact Finset.mem_univ v)

@[simp] theorem archADiv_apply {X : Scheme.{0}} (g : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) (v : InfinitePlace F) :
    archADiv F g xF v = (InfinitePlace.mult v : ℝ) * g (archPoint xF v) := rfl

/-- ★★★**Green 関数の引き戻しは加法的**——**仮定を 1 つも要しない**。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★★テンソル積で計量が掛かる(= 対数が足される)ことの、引き戻し側での帰結である。
★有限素点側と違って、ここには `comap` の穴が無い。 -/
theorem archADiv_add {X : Scheme.{0}} (g h : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) :
    archADiv F (fun p => g p + h p) xF = archADiv F g xF + archADiv F h xF := by
  ext v
  simp only [archADiv_apply, Finsupp.add_apply]
  ring

/-- ★**零 Green 関数の引き戻しは零**。 -/
@[simp] theorem archADiv_zero {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    archADiv F (fun _ : complexPoints X => (0 : ℝ)) xF = 0 := by
  ext v
  simp only [archADiv_apply, Finsupp.coe_zero, Pi.zero_apply, mul_zero]

/-! ## ★★★算術因子の引き戻し —— 2 成分を束ねる -/

/-- ★★**算術因子**(有効 Cartier 因子 + Green 関数)。

★原文の算術直線束 `L̄ = (L, |−|_L)` を**因子の側で**表したものである。
`Pic(X) ≅ CaCl(X)` を使って可逆層を避けている——
mathlib に可逆層も `SheafOfModules` のテンソル積も無いからである(2026-08-16 実測)。 -/
structure ArithCartier (X : Scheme.{0}) where
  /-- 有限素点側: 有効 Cartier 因子(イデアル層) -/
  divisor : X.IdealSheafData
  /-- 無限素点側: Green 関数 -/
  green : GreenFn X

/-- ★★**テンソル積** —— 因子は掛け、Green 関数は足す。

★★この非対称(積と和)が、対数を取って持つ理由である。 -/
def ArithCartier.tensor {X : Scheme.{0}} (D E : ArithCartier X) : ArithCartier X :=
  { divisor := D.divisor * E.divisor
    green := fun p => D.green p + E.green p }

/-- ★**空の算術因子**(単位元)。 -/
def ArithCartier.one (X : Scheme.{0}) : ArithCartier X :=
  { divisor := ⊤, green := fun _ => 0 }

/-- ★★**算術因子の引き戻し** —— `x_F^* D̄ : ADiv F`。

有限素点側は `pullbackIdeal` → `idealADiv`、無限素点側は `archADiv`。 -/
noncomputable def pullbackADiv {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) : ADiv F :=
  ((idealADiv F (pullbackIdeal F D.divisor xF)).fin, archADiv F D.green xF)

@[simp] theorem pullbackADiv_fin {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) :
    (pullbackADiv F D xF).fin = (idealADiv F (pullbackIdeal F D.divisor xF)).fin := rfl

@[simp] theorem pullbackADiv_arc {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) :
    (pullbackADiv F D xF).arc = archADiv F D.green xF := rfl

/-- ★**単位イデアルの算術因子は零**。

★★証明は `idealADiv_mul` の**再利用**である——
`⊤ · ⊤ = ⊤` から `a = a + a`、したがって `a = 0`。
分解を直接計算する必要が無い。 -/
theorem idealADiv_top : idealADiv F (⊤ : Ideal (𝓞 F)) = 0 := by
  have hne : (⊤ : Ideal (𝓞 F)) ≠ 0 := by
    intro h
    exact (one_ne_zero (α := Ideal (𝓞 F))) (by rw [Ideal.one_eq_top]; exact h)
  have h := idealADiv_mul F ⊤ ⊤ hne hne
  rw [Ideal.mul_top] at h
  have h0 : idealADiv F (⊤ : Ideal (𝓞 F)) + 0
      = idealADiv F (⊤ : Ideal (𝓞 F)) + idealADiv F (⊤ : Ideal (𝓞 F)) := by
    rw [add_zero]; exact h
  exact (add_left_cancel h0).symm

/-- ★**空の算術因子の引き戻しは零**(非空虚性の witness)。 -/
@[simp] theorem pullbackADiv_one {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    pullbackADiv F (ArithCartier.one X) xF = 0 := by
  refine Prod.ext ?_ ?_
  · show (pullbackADiv F (ArithCartier.one X) xF).fin = (0 : ADiv F).fin
    rw [pullbackADiv_fin]
    show (idealADiv F (pullbackIdeal F (⊤ : X.IdealSheafData) xF)).fin = _
    rw [pullbackIdeal_top, idealADiv_top]
  · show (pullbackADiv F (ArithCartier.one X) xF).arc = (0 : ADiv F).arc
    rw [pullbackADiv_arc]
    exact archADiv_zero F xF

/-- ★★★**引き戻しの無限素点側はテンソル積で加法的** —— 無条件。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★★`Proposition 1.4, (i)` の**アルキメデス半分**である。
★有限素点半分は `IdealSheafData.comap` の積の保存を要し、そこは穴のままである。 -/
theorem pullbackADiv_arc_tensor {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) :
    (pullbackADiv F (D.tensor E) xF).arc
      = (pullbackADiv F D xF).arc + (pullbackADiv F E xF).arc := by
  rw [pullbackADiv_arc, pullbackADiv_arc, pullbackADiv_arc]
  exact archADiv_add F D.green E.green xF

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.1` 全体には可逆層のテンソル積が要り、
`Proposition 1.4, (i)` 全体には有限素点側の積の保存が要る。 -/

def archPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(ℂ-点と Green 関数の引き戻しのみ)",
    sectionId := "genell-def-1-1-i" }

def pullbackADiv_arc_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Proposition 1.4, (i)(アルキメデス側のみ——Green 関数の加法性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.GenEll
