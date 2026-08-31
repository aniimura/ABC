/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AmpleImmersion
import ABC3.Found.GenEll.HyperplanePullbackGlobal
import ABC3.Found.GenEll.NorthcottVeryAmple
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★高さの側 —— 貼った射に沿った Northcott（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6-7。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★これは何か —— 幾何と高さが出会う

`§9-920` で `ψ : X ⟶ ℙᴺ` が `IsAmple` から作れた。★本ファイルは**高さの側**を繋ぐ:

| 段 | 内容 |
|---|---|
| 1 | `ψ` に沿った超平面の引き戻しは `(s_0/s_i)(x)` が生成する |
| 2 | `ht_{ψ^*Ē}(x) = ht_Ē(ψ ∘ x)`（関手性、`§9-879`） |
| 3 | ★**`ht_{ψ^*Ē}(x) = log H(座標)/[F:ℚ]`**——素朴高さそのもの |
| 4 | ★★`northcott_comap`（`§9-881`）の `hcomp` が**チャートを通る点については自動で埋まる** |

★★★これで `northcott_of_veryAmple`（`§9-882`）の仮定のうち
**`hcomp`（点がチャートを通ること）が `ψ` の構成から出る**ようになった。

## ★★機構 —— 貼った射はチャートで元の射に戻る

`ι_globalToProj`（`§9-911`）が

    `X_{s_i} ↪ X ⟶ ℙᴺ  =  globalChartToProj i`

を与えるので、`§9-886`（チャート射に沿った超平面の引き戻し）と
`§9-871`（超平面の算術因子の高さは素朴高さ）が**そのまま当たる**。

## ★残っている段（明示）

★★★★残るのは **`D^{⊗n} = ψ^*Ē`**（与えられた算術因子と引き戻しの一致）であり、
これは段 E0（切断の零点因子 `divisorOfSection`）の内容である。
★`Scheme.IdealSheafData.ext_of_iSup_eq_top`（mathlib）が
「アフィン被覆の上で一致すれば等しい」を与えるので、
`§9-885`（チャートの上で超平面の引き戻しは `(s_0/s_i)`）を貼れば出る見込みである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★段 1 —— 引き戻しイデアル -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**`ψ` に沿った超平面の引き戻しは `(s_0/s_i)(x)` が生成する**。

★`ι_globalToProj`（`§9-911`）で `§9-886` に戻すだけである。 -/
theorem pullbackIdeal_hyperplane_globalToProj (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme) :
    pullbackIdeal F (hyperplaneIdeal N ℤ)
        (xF ≫ ((nonVanishing M (s i)).ι ≫ globalToProj M hM φ s hcov))
      = Ideal.span {(Spec.preimage (xF ≫ (nonVanishing M (s i)).toScheme.toSpecΓ)).hom
          (((nonVanishing M (s i)).topIso.inv.hom) (globalRatio M hM (s 0) (s i)))} := by
  rw [ι_globalToProj]
  exact pullbackIdeal_hyperplane_globalRatio M hM φ s i F xF

/-! ## ★★★★★★★★★★★★段 2-3 —— 高さは素朴高さである -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★**関手性をチャートで読む**。 -/
theorem htArith_hyperplaneArith_comap (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme) :
    htArith F ((hyperplaneArith N).comap (globalToProj M hM φ s hcov))
        (xF ≫ (nonVanishing M (s i)).ι)
      = htArith F (hyperplaneArith N) (xF ≫ globalChartToProj M hM φ s i) := by
  rw [htArith_comap, Category.assoc, ι_globalToProj]

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★**`ht_{ψ^*Ē}(x)` は座標の素朴高さである**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★これが幾何（`§9-920` の `ψ`）と高さ（`§9-871`）の**出会う点**である。 -/
theorem htArith_globalToProj_eq_log_mulHeight (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme)
    (x : Fin (N + 1) → NumberField.RingOfIntegers F) (hx : x ≠ 0)
    (hw : ∀ k, x k = (Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)).hom
      (projCoord N ℤ i k) * x i)
    (h0 : x 0 ≠ 0) :
    htArith F ((hyperplaneArith N).comap (globalToProj M hM φ s hcov))
        (xF ≫ (nonVanishing M (s i)).ι)
      = Real.log (Height.mulHeight (fun k => ((x k : F)))) / (Module.finrank ℚ F : ℝ) := by
  rw [htArith_hyperplaneArith_comap]
  have hfac : xF ≫ globalChartToProj M hM φ s i
      = Spec.map (Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)) ≫ chartA N ℤ i := by
    rw [Spec.map_preimage, globalChartToProj, Category.assoc]
  rw [hfac]
  exact htArith_hyperplaneArith F N i _ x hx hw h0

/-! ## ★★★★★★★★★★★★★★★★段 4 —— `hcomp` が自動で埋まる -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★**チャートを通る点は `Spec.map ρ ≫ chartA` の形に分解する**。 -/
theorem chart_factorization (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme) :
    (xF ≫ (nonVanishing M (s i)).ι) ≫ globalToProj M hM φ s hcov
      = Spec.map (Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)) ≫ chartA N ℤ i := by
  rw [Category.assoc, ι_globalToProj, Spec.map_preimage, globalChartToProj, Category.assoc]

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★**貼った射に沿った Northcott**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★`northcott_comap`（`§9-881`）の `hcomp`（点がチャートを通ること）が
**`ψ` の構成から自動で埋まる**——それが `chart_factorization` である。

★★★残る仮定は
* 点がどれかのチャートを通ること（`chart`・`xF`）——`X` が固有なら自動のはず
* 座標が点を分けること（`hinj`）——`ψ` が埋め込み（`§9-920`）であることの点版 -/
theorem northcott_globalToProj (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ} (d : ℕ)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p))
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (chart : P → Fin (N + 1))
    (xF : ∀ p, haveI := hnf p;
      specRingOfIntegers (fld p) ⟶ (nonVanishing M (s (chart p))).toScheme)
    (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
    (hx : ∀ p, x p ≠ 0)
    (hw : ∀ p, haveI := hnf p; ∀ k, x p k =
      (Spec.preimage (xF p ≫ globalChartMorphism M hM φ s (chart p))).hom
        (projCoord N ℤ (chart p) k) * x p (chart p))
    (h0 : ∀ p, x p 0 ≠ 0)
    (idx : Fin (N + 1))
    (hinj : Function.Injective (fun (p : P) (k : Fin (N + 1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
    (C : ℝ) :
    {p : P | haveI := hnf p;
      htArith (fld p) (E.comap (globalToProj M hM φ s hcov))
        (xF p ≫ (nonVanishing M (s (chart p))).ι) ≤ C}.Finite := by
  refine northcott_comap N d E hdiv hcont (globalToProj M hM φ s hcov) fld hnf hdeg
    (fun p => xF p ≫ (nonVanishing M (s (chart p))).ι) x (fun p => ?_) idx hinj C
  haveI := hnf p
  exact ⟨chart p, Spec.preimage (xF p ≫ globalChartMorphism M hM φ s (chart p)),
    chart_factorization M hM φ s hcov (chart p) (fld p) (xF p), hx p, hw p, h0 p⟩

/-! ## ★出典の紐付け(`.src`) -/

def pullbackIdeal_hyperplane_globalToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(ψ に沿った超平面の引き戻しは (s_0/s_i)(x) が生成する)",
    sectionId := "genell-prop-1-4" }

def htArith_hyperplaneArith_comap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(高さの関手性をチャートで読む)",
    sectionId := "genell-prop-1-4" }

def htArith_globalToProj_eq_log_mulHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ψ^*Ē の高さは座標の素朴高さである)",
    sectionId := "genell-prop-1-4" }

def chart_factorization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートを通る点の分解)",
    sectionId := "genell-prop-1-4" }

def northcott_globalToProj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(貼った射に沿った Northcott)",
    sectionId := "genell-prop-1-4" }

def northcott_globalToProj.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_comap(射に沿った Northcott の移送、§9-881)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_comap") 3,
    .citation "[ABC3]" "ι_globalToProj(貼った射はチャートで元の射に戻る、§9-911)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ι_globalToProj") 2,
    .citation "[ABC3]" "htArith_hyperplaneArith(超平面の算術因子の高さは素朴高さ、§9-871)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_hyperplaneArith") 3,
    .implicitStep
      ("★★★これで northcott_of_veryAmple(§9-882)の仮定のうち " ++
       "**hcomp(点がチャートを通ること)が ψ の構成から出る**ようになった") 3,
    .implicitStep
      ("★★★★残るのは D^{⊗n} = ψ^*Ē(与えられた算術因子と引き戻しの一致)であり、" ++
       "これは段 E0(切断の零点因子 divisorOfSection)の内容である。" ++
       "★Scheme.IdealSheafData.ext_of_iSup_eq_top(mathlib)が" ++
       "「アフィン被覆の上で一致すれば等しい」を与えるので、" ++
       "§9-885(チャートの上で超平面の引き戻しは (s_0/s_i))を貼れば出る見込みである") 5 ]

end ABC3.Found.GenEll
