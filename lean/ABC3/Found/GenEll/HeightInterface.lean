import ABC3.Found.GenEll.UPoint
import ABC3.Interface.GenEll.ArithLineBundle

/-!
# [GenEll] Definition 1.1, (i) —— **`PulledBackClassData` の非空虚 witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3–p.4。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★`waiting` を外す

`Interface/GenEll/ArithLineBundle.lean` は `PulledBackClassData` に
**`waiting` を置いていた**——「退化 witness(すべてを `0` に送る操作)なら今すぐ作れるが、
それで G2 を満たしても作業キューから消えるだけで何も進まない」という理由である。

★★★**本ファイルはその退化 witness ではない。**供給するのは

| 場 | 中身 |
|---|---|
| `Point` | `UPoint X D`(底変換で同一視した代数的点。`UPoint.lean`) |
| `Bundle` | `{g : GreenFn X // IsConjInvariant g}`(原文の「`ι_X` と両立する計量」) |
| `degOver F` | **`F` 上のモデルを 1 つ選んで `deg_F(x_F^* M̄)` を計算したもの** |
| `height` | `htU`(`Quot.lift` で商に降ろした高さ) |
| `base_change_invariant` | ★**定理**(`htArith_baseChange_natural` から出る) |

★★**要である `base_change_invariant` を仮説ではなく定理として与えている**——
そこが退化 witness との違いである。原文 p.4 の
`deg_K(L̄|_{Spec 𝒪_K}) = deg_F(L̄)` を、我々は**証明した**。

## ★★★何を供給していないか(条を明示する)

- `Point` は `U_X(ℚ̄) = X(ℚ̄) \ D` であって **`X(ℚ̄)` 全体ではない**
  (因子表示なので「`x` が `D` を通らない」が要る)。
- `Bundle` は**因子 `D` を固定し計量だけを動かした族**であって、
  一般の可逆層ではない。

★★★**それでも `PulledBackClassData` の witness としては十分である**——
この `structure` は設計上「算術直線束そのもの」ではなく
**高さの定義が実際に使う操作だけ**を受けているからである
(`Interface` 側の docstring がそう書いている)。

★点の型が**空でない**ことは `Check/GenEll/HeightInterfaceNondegenerate.lean` で確かめる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField
open ABC3.Interface.GenEll

variable {X : Scheme.{0}}

/-! ## ★`x` が `F` 上で定義されること -/

/-- ★**原文の `x ∈ X(F)`** —— `F` 上のモデルが存在すること。 -/
def DefinedOverField (D : X.IdealSheafData) (F : Type) [Field F] [NumberField F]
    (x : UPoint X D) : Prop :=
  ∃ (xF : specRingOfIntegers F ⟶ X) (hJ : pullbackIdeal F D xF ≠ 0),
    UPoint.mk (algPointOff F xF hJ) = x

/-- ★`F` 上のモデルからは、それが `F` 上で定義されることが出る。 -/
theorem definedOverField_mk (D : X.IdealSheafData) (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X) (hJ : pullbackIdeal F D xF ≠ 0) :
    DefinedOverField D F (UPoint.mk (algPointOff F xF hJ)) :=
  ⟨xF, hJ, rfl⟩

/-! ## ★★`deg_F(x_F^* M̄)` -/

open scoped Classical in
/-- ★★**原文 p.4 の `deg_F(x_F^* M̄)`** —— `F` 上のモデルを 1 つ選んで計算する。

★★★**選び方に依らないことが `degOverField_eq` である**——
そこが内容であって、`Classical.choose` は飾りである。 -/
noncomputable def degOverField (D : ArithCartier X) (F : Type) [Field F] [NumberField F]
    (x : UPoint X D.divisor) : ℝ :=
  if h : DefinedOverField D.divisor F x then htArith F D h.choose else 0

/-- ★★★**選んだモデルに依らない**——値は `htU` に一致する。

★これが `Interface` の `height_eq` そのものであり、
**`base_change_invariant` もここから出る**(2 つの `F`・`K` で同じ値を経由する)。 -/
theorem degOverField_eq (D : ArithCartier X) (hg : IsConjInvariant D.green)
    (F : Type) [Field F] [NumberField F] (x : UPoint X D.divisor)
    (h : DefinedOverField D.divisor F x) :
    degOverField D F x = htU D hg x := by
  obtain ⟨hJ, hmk⟩ := h.choose_spec
  calc degOverField D F x = htArith F D h.choose := by rw [degOverField, dif_pos h]
    _ = htOff D (algPointOff F h.choose hJ) := rfl
    _ = htU D hg (UPoint.mk (algPointOff F h.choose hJ)) := (htU_mk D hg _).symm
    _ = htU D hg x := by rw [hmk]

/-! ## ★★★`Interface` への供給 -/

/-- ★★★**`PulledBackClassData` の実例**——因子 `D₀` を固定し、計量を動かした族。

★`degOver` は `F` 上のモデルを選んで計算した `deg_F(x_F^* M̄)` であり、
**`0` に潰した退化 witness ではない**。
★★`base_change_invariant` は `degOverField_eq` を 2 度使って出る——
**仮説ではなく定理である**。 -/
noncomputable def heightData (D₀ : X.IdealSheafData) :
    PulledBackClassData.{1, 0} where
  Point := UPoint X D₀
  Bundle := {g : GreenFn X // IsConjInvariant g}
  DefinedOver := fun F _ _ x => DefinedOverField D₀ F x
  degOver := fun F _ _ g x => degOverField ⟨D₀, g.1⟩ F x
  height := fun g x => htU ⟨D₀, g.1⟩ g.2 x
  base_change_invariant := fun F K _ _ _ _ _ g x hF hK => by
    rw [degOverField_eq ⟨D₀, g.1⟩ g.2 K x hK, degOverField_eq ⟨D₀, g.1⟩ g.2 F x hF]
  height_eq := fun F _ _ g x h => (degOverField_eq ⟨D₀, g.1⟩ g.2 F x h).symm

/-! ## ★★点の型が空でないこと -/

/-- ★★**`Spec 𝓞_ℚ` 自身の上の空因子には点がある**——恒等射である。

★★★これがあるので `heightData` は**空虚な witness ではない**
(`Point` が空なら `∀` はすべて自明に成り立ってしまう)。 -/
theorem nonempty_uPoint_specRingOfIntegers :
    Nonempty (UPoint (specRingOfIntegers ℚ) (⊤ : (specRingOfIntegers ℚ).IdealSheafData)) := by
  refine ⟨UPoint.mk (algPointOff ℚ (𝟙 _) ?_)⟩
  rw [pullbackIdeal_top]
  exact fun h => absurd (h ▸ Submodule.mem_top (x := (1 : 𝓞 ℚ))) (by simp)

end ABC3.Found.GenEll

namespace ABC3.Interface.GenEll

open AlgebraicGeometry NumberField

/-- ★★★**G2 非空虚 witness**。`sorry` 無し。

`structure` は `Interface/GenEll/ArithLineBundle.lean` にあるが、
witness は実装を要するのでこちら(`Found/`)で `Interface.GenEll` 名前空間に足す
(`Found/PGC/ResidueCardinality.lean` と同じ向き)。

★★★**`waiting` が言っていた「層 B(可逆層)・層 C(解析化)」は
まだ無い。**しかし `PulledBackClassData` は設計上それらを型に出しておらず、
**高さの定義が実際に使う操作だけ**を受けている。その操作は
`Found/GenEll/` で構成済みである(因子表示・`U_X(ℚ̄)` の範囲)。

★点の型が空でないことは
`ABC3.Found.GenEll.nonempty_uPoint_specRingOfIntegers`。 -/
theorem PulledBackClassData.nonvacuous :
    Nonempty PulledBackClassData.{1, 0} :=
  ⟨ABC3.Found.GenEll.heightData
    (X := ABC3.Found.GenEll.specRingOfIntegers ℚ) ⊤⟩

end ABC3.Interface.GenEll

namespace ABC3.Found.GenEll

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文の `Definition 1.1, (i)` は**可逆層と hermitian 計量**で
算術直線束を定めるが、本ファイルが供給するのは**因子と Green 関数**の表示である。 -/

def heightData.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(高さの定義が使う操作のみ——因子表示・U_X(ℚ̄) の範囲)",
    sectionId := "genell-def-1-1-i" }

def degOverField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(deg_F(x_F^* M̄) がモデルの取り方に依らないこと)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
