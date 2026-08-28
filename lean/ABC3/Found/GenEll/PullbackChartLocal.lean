/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottGlobalToProj
import ABC3.Found.GenEll.PullbackLocalization
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★引き戻しイデアルを**素点ごとに**読む（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★これは何か —— `chart` 仮定の**幾何の側**の第 1 歩

`§9-928` で測ったとおり、点 `Spec 𝓞_F ⟶ ℙᴺ_ℤ` が**大域的に**チャート `D₊(x_i)` を
通るのは「座標の生成するイデアル `𝔞` が `x_i` で生成される」ときに限る
——**イデアル類が自明でないと通らない**。

★★しかし**局所的には必ず通る**: 素イデアル `𝔭` ごとに `v_𝔭(x_i)` が最小になる `i`
を取れば `𝔞_𝔭 = (x_i)_𝔭` であり、`Spec (𝓞_F)_𝔭` はチャートを通る。

★★★★**本ファイルはそのための配線**である:

1. `pullbackIdealOf` は**任意の底環 `B`** で動く（`§9-865`・`§9-907` は `B = 𝓞_F` だった）
2. `pullbackIdealOf B (D.comap f) x = pullbackIdealOf B D (x ≫ f)`（合成）
3. `Spec B ⟶ X_{s_i}` を通る点では引き戻しは `(s_0/s_i)(x)` が生成する
4. ★**局所化した点がチャートを通れば、引き戻しイデアルの局所化が読める**

## ★★★機構 —— 既にある 2 本を繋ぐだけである

* `pullbackIdealOf_specMap`（`§9-907`）——任意の環準同型に沿って拡大イデアルになる
* `pullbackIdeal_map_algebraMap`（同）——`𝓞_F → A` への拡大は `Spec A` の点の引き戻し

★★どちらも**数体を仮定していない**ので、`§9-886`・`§9-921` の証明の骨格が
`B = (𝓞_F)_𝔭` でもそのまま通る。

## ★残っている段（明示）

★★★残るのは**算術の側**——局所生成元 `g_𝔭` が `g_𝔭·x_{i(𝔭)} = x_0` を満たすことから

    `pullbackIdeal · 𝔞 = (x_0)`

を出す段である（`Ideal.eq_of_localization_maximal`）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★合成に沿った引き戻し -/

/-- ★★★**`comap` に沿った引き戻しは合成の引き戻しである**。

★`Scheme.IdealSheafData.comap_comp` をそのまま `pullbackIdealOf` の定義に流すだけ。 -/
theorem pullbackIdealOf_comap {B : CommRingCat.{0}} {X Y : Scheme.{0}}
    (D : Y.IdealSheafData) (f : X ⟶ Y) (x : Spec B ⟶ X) :
    pullbackIdealOf B (D.comap f) x = pullbackIdealOf B D (x ≫ f) := by
  rw [pullbackIdealOf, pullbackIdealOf, Scheme.IdealSheafData.comap_comp]

/-! ## ★★★★★★★★★★★★任意の底環の上でのチャートの読み -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★**チャート射に沿った超平面の引き戻し**——任意の底環 `B` で。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-886` は `B = 𝓞_F` の場合であった。証明は同じ
——`globalChartMorphism` の終域がアフィンなので `Spec.map_preimage` が効く。 -/
theorem pullbackIdealOf_hyperplane_globalChart (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (B : CommRingCat.{0})
    (xB : Spec B ⟶ (nonVanishing M (s i)).toScheme) :
    pullbackIdealOf B (hyperplaneIdeal N ℤ) (xB ≫ globalChartToProj M hM φ s i)
      = Ideal.span {(Spec.preimage (xB ≫ globalChartMorphism M hM φ s i)).hom
          (projCoord N ℤ i 0)} := by
  have h : xB ≫ globalChartToProj M hM φ s i
      = Spec.map (Spec.preimage (xB ≫ globalChartMorphism M hM φ s i)) ≫ chartA N ℤ i := by
    rw [Spec.map_preimage, globalChartToProj, Category.assoc]
  rw [h, pullbackIdealOf_hyperplane_point]

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★**生成元は比の切断の値である**——任意の底環 `B` で。

★`§9-886` の `preimage_globalChartMorphism_projCoord` の一般化。 -/
theorem preimage_globalChartMorphism_projCoord' (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1))
    (B : CommRingCat.{0})
    (xB : Spec B ⟶ (nonVanishing M (s i)).toScheme) :
    (Spec.preimage (xB ≫ globalChartMorphism M hM φ s i)).hom (projCoord N ℤ i j)
      = (Spec.preimage (xB ≫ (nonVanishing M (s i)).toScheme.toSpecΓ)).hom
          (((nonVanishing M (s i)).topIso.inv.hom) (globalRatio M hM (s j) (s i))) := by
  rw [globalChartMorphism, ← Category.assoc, specPreimage_comp_specMap]
  show (Spec.preimage (xB ≫ (nonVanishing M (s i)).toScheme.toSpecΓ)).hom
    (((nonVanishing M (s i)).topIso.inv.hom)
      (globalAwayHom M hM φ s i (projCoord N ℤ i j))) = _
  rw [globalAwayHom_projCoord]

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★**`ψ` に沿った超平面の引き戻し**——任意の底環 `B` で。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-921` の `pullbackIdeal_hyperplane_globalToProj` の一般化である。 -/
theorem pullbackIdealOf_hyperplane_globalToProj (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1))
    (B : CommRingCat.{0})
    (xB : Spec B ⟶ (nonVanishing M (s i)).toScheme) :
    pullbackIdealOf B (hyperplaneIdeal N ℤ)
        (xB ≫ ((nonVanishing M (s i)).ι ≫ globalToProj M hM φ s hcov))
      = Ideal.span {(Spec.preimage (xB ≫ (nonVanishing M (s i)).toScheme.toSpecΓ)).hom
          (((nonVanishing M (s i)).topIso.inv.hom) (globalRatio M hM (s 0) (s i)))} := by
  rw [ι_globalToProj, pullbackIdealOf_hyperplane_globalChart,
    preimage_globalChartMorphism_projCoord']

/-! ## ★★★★★★★★★★★★★★★局所化した点で読む -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★**局所化した点がチャートを通れば、引き戻しイデアルの
局所化が読める**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `Spec A ⟶ Spec 𝓞_F ⟶ X` が `X_{s_i}` を通る
      ⟹ `(ψ^*超平面 の引き戻し)·A = ( (s_0/s_i)(y) )`

★★★これが `chart` 仮定を**大域から局所へ**降ろす段である
——`A = (𝓞_F)_𝔭` を取れば、素点ごとには必ずチャートが取れる。 -/
theorem map_pullbackIdeal_globalToProj_of_localChart (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (A : Type) [CommRing A] [Algebra (NumberField.RingOfIntegers F) A]
    (xF : specRingOfIntegers F ⟶ X)
    (y : Spec (CommRingCat.of A) ⟶ (nonVanishing M (s i)).toScheme)
    (hy : Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF
        = y ≫ (nonVanishing M (s i)).ι) :
    (pullbackIdeal F ((hyperplaneIdeal N ℤ).comap (globalToProj M hM φ s hcov)) xF).map
        (algebraMap (NumberField.RingOfIntegers F) A)
      = Ideal.span {(Spec.preimage (y ≫ (nonVanishing M (s i)).toScheme.toSpecΓ)).hom
          (((nonVanishing M (s i)).topIso.inv.hom) (globalRatio M hM (s 0) (s i)))} := by
  rw [pullbackIdeal_map_algebraMap, pullbackIdealOf_comap, hy, Category.assoc,
    pullbackIdealOf_hyperplane_globalToProj]

/-! ## ★出典の紐付け(`.src`) -/

def pullbackIdealOf_comap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(comap に沿った引き戻しは合成の引き戻しである)",
    sectionId := "genell-def-1-2-i" }

def pullbackIdealOf_hyperplane_globalChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャート射に沿った超平面の引き戻し——任意の底環で)",
    sectionId := "genell-prop-1-4" }

def preimage_globalChartMorphism_projCoord'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(生成元は比の切断の値である——任意の底環で)",
    sectionId := "genell-prop-1-4" }

def pullbackIdealOf_hyperplane_globalToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ψ に沿った超平面の引き戻し——任意の底環で)",
    sectionId := "genell-prop-1-4" }

def map_pullbackIdeal_globalToProj_of_localChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(局所化した点がチャートを通れば引き戻しイデアルの局所化が読める)",
    sectionId := "genell-prop-1-4" }

def map_pullbackIdeal_globalToProj_of_localChart.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullbackIdeal_map_algebraMap(拡大イデアルは Spec A の点の引き戻し、§9-907)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackIdeal_map_algebraMap") 2,
    .citation "[ABC3]" "pullbackIdealOf_specMap(任意の環準同型に沿って拡大イデアルになる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackIdealOf_specMap") 2,
    .citation "[ABC3]" "ι_globalToProj(貼った射のチャート上での値、§9-911)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ι_globalToProj") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): 点 Spec 𝓞_F ⟶ ℙᴺ が**大域的に**チャート D₊(x_i) を通るのは" ++
       "イデアル類が自明なときに限る(§9-928)が、**局所的には必ず通る**" ++
       "——素イデアル 𝔭 ごとに v_𝔭(x_i) が最小になる i を取れば 𝔞_𝔭 = (x_i)_𝔭 である") 5,
    .implicitStep
      ("★★§9-865・§9-886・§9-921 の証明は**数体を仮定していない**ので、" ++
       "底環を B = (𝓞_F)_𝔭 に替えてもそのまま通る" ++
       "——pullbackIdealOf_specMap が一般の可換環で成り立つからである") 3,
    .implicitStep
      ("★★★残るのは算術の側——局所生成元 g_𝔭 が g_𝔭·x_{i(𝔭)} = x_0 を満たすことから " ++
       "pullbackIdeal · 𝔞 = (x_0) を出す段である(Ideal.eq_of_localization_maximal)") 3 ]

end ABC3.Found.GenEll
