import ABC3.Found.Arakelov.ArcEvalOn

/-!
# Arakelov (C3) 第 257 ブロック —— ★★★★★★**生成切断で正規化したノルム**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★選択を**比で消す**

第 248 の `arcMetricOf` は各点で自明化を**勝手に選ぶ**ので、
ノルム自体は `p` について暴れる(連続でない)。

★★★しかし**比**を取れば選択は消える:

    ‖w‖ := nrm(w) / nrm(g(p))          (g は `V` 上の生成切断)

★1 次元ファイバー上のノルムは正の定数倍しか違わないので、
分子と分母で**同じ定数が掛かって約分される**。

★★★★★したがって `genNorm` は**選択に依らず**、`s = c·g` のとき

    genNorm (arcEvalOnTop p V h s) = |c(p)|

となる——**連続性が正則関数の連続性に落ちる**。

## ★★これが §9-297 の壁を越えた後の主要部品である

| 定義・定理 | 内容 |
|---|---|
| `arcEvalOn_smul` | ★評価は係数環上線形(押し出しの綴りで述べる) |
| `refNorm` | ★基準ノルム(選択に依る) |
| `genNorm` | ★★★★**正規化したノルム** |
| `genNorm_nonneg` / `_smul` / `_eq_zero_iff` / `_self` | ★★4 法則 |

★`genNorm_self : genNorm (g(p)) = 1` が正規化の意味である。

## ★摩擦 —— 線形性は**押し出しの綴り**で述べる

`arcEvalOn` の終域を引き戻し側 `((pullback p).obj F).val.obj (op (p⁻¹U))` で書くと
`HSMul Γ(X,U) …` が解決しない。★**押し出し側** `(p_* p^* F).val.obj (op U)` で
述べると解決する——加群構造がそちらに登録されているからである。
★★[[ring-instance-two-paths]]の「**登録されている側で述べる**」型の逃げ道。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}} (F : X.Modules) (p : Spec (CommRingCat.of ℂ) ⟶ X)

/-- ★`arcEvalOn` は係数環上線形である(押し出しの綴りで述べる)。 -/
theorem arcEvalOn_smul (U : X.Opens) (c : ((X.presheaf.obj (op U)) : Type))
    (s : (F.val.obj (op U) : Type)) :
    (((((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app F).val.app (op U)).hom (c • s)) :
        ((((Scheme.Modules.pushforward p).obj
          ((Scheme.Modules.pullback p).obj F)).val.obj (op U)) : Type))
      = c • ((((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app F).val.app (op U)).hom s) :=
  map_smul _ c s


variable (hF : IsLocallyTrivial X F.val)

/-- ★★基準ノルム(選択に依る)。 -/
noncomputable def refNorm (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) : ℝ :=
  (arcMetricOf F hF).nrm p w

/-- ★★★★生成切断で正規化したノルム——**選択に依らない**。 -/
noncomputable def genNorm (V : X.Opens) (g : (F.val.obj (op V) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤) (w : ↥(arcFiber p F)) : ℝ :=
  refNorm F hF p w / refNorm F hF p (arcEvalOnTop F p V h g)

theorem genNorm_nonneg (V : X.Opens) (g : (F.val.obj (op V) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤) (w : ↥(arcFiber p F)) :
    0 ≤ genNorm F hF V g p h w :=
  div_nonneg ((arcMetricOf F hF).nonneg p w) ((arcMetricOf F hF).nonneg p _)

theorem genNorm_smul (V : X.Opens) (g : (F.val.obj (op V) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤)
    (c : (CommRingCat.of ℂ : Type)) (w : ↥(arcFiber p F)) :
    genNorm F hF V g p h (c • w) = ‖c‖ * genNorm F hF V g p h w := by
  show (arcMetricOf F hF).nrm p (c • w) / _ = _
  rw [(arcMetricOf F hF).smul p c w]
  exact mul_div_assoc _ _ _

theorem genNorm_eq_zero_iff (V : X.Opens) (g : (F.val.obj (op V) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤)
    (hg : arcEvalOnTop F p V h g ≠ 0) (w : ↥(arcFiber p F)) :
    genNorm F hF V g p h w = 0 ↔ w = 0 := by
  have hd : refNorm F hF p (arcEvalOnTop F p V h g) ≠ 0 := by
    intro hz
    exact hg (((arcMetricOf F hF).eq_zero_iff p _).1 hz)
  show (arcMetricOf F hF).nrm p w / refNorm F hF p (arcEvalOnTop F p V h g) = 0 ↔ w = 0
  rw [div_eq_zero_iff]
  constructor
  · rintro (hw | hz)
    · exact ((arcMetricOf F hF).eq_zero_iff p w).1 hw
    · exact absurd hz hd
  · intro hw
    exact Or.inl (((arcMetricOf F hF).eq_zero_iff p w).2 hw)


theorem genNorm_self (V : X.Opens) (g : (F.val.obj (op V) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (h : p ⁻¹ᵁ V = ⊤)
    (hg : arcEvalOnTop F p V h g ≠ 0) :
    genNorm F hF V g p h (arcEvalOnTop F p V h g) = 1 := by
  have hd : refNorm F hF p (arcEvalOnTop F p V h g) ≠ 0 := by
    intro hz
    exact hg (((arcMetricOf F hF).eq_zero_iff p _).1 hz)
  exact div_self hd


/-! ## ★出典の紐付け(`.src`) -/

def genNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——生成切断で正規化したノルム)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
