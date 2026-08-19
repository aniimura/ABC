import ABC3.Found.Arakelov.ArcEvalGlobal

/-!
# Arakelov (C3) 第 253 ブロック —— ★★★★★★**自明な直線束には連続な計量が入る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★第 244–252 の合流点

自明化 `t : L ≅ 𝒪_V` があれば、ファイバーの同型

    arcFiber p L ≅ arcFiber p 𝒪_V ≅ Γ(Spec ℂ, ⊤) ≅ ℂ

を通してノルムが入る。★★そして**切断のノルムは正則関数の絶対値**である:

    ‖s‖(p) = |(t s)(p)|            (`trivNorm_arcEval`)

★★★右辺は第 252(大域の正則関数は連続)でそのまま連続——**連続性が閉じた**。

## ★★これで「自明な場合」の計量は完成である

| 法則 | 定理 |
|---|---|
| 非負 | `trivNorm_nonneg` |
| `0 ↔ v = 0` | `trivNorm_eq_zero_iff` |
| `‖c·v‖ = \|c\|‖v‖` | `trivNorm_smul` |
| ★**連続性** | `continuous_trivNorm` |
| まとめ | ★★★`trivArcMetric` |

★★★★残るのは**貼り合わせ**(1 の分割)だけである。

## ★摩擦 —— `rw` を 3 回やめた

`presheaf` の二重路(`TopCat.Presheaf` vs `_ᵒᵖ ⥤ _`)で `rw [← hnat]` が落ちた。
★★逃げ道は毎回同じ:**`congrArg` で外側を包んで `trans` で繋ぐ**。
★★★`simpa` も「簡約後の型が合わない」で落ちるので、`calc` で書き下した。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {V : Scheme.{0}} {L : V.Modules}

/-- ★自明化からファイバーの同型を作る。 -/
noncomputable def trivFiberIso (t : L ≅ unitModules V) (p : Spec (CommRingCat.of ℂ) ⟶ V) :
    arcFiber p L ≅ moduleSpecΓFunctor (R := CommRingCat.of ℂ) |>.obj
      (unitModules (Spec (CommRingCat.of ℂ))) :=
  (moduleSpecΓFunctor (R := CommRingCat.of ℂ)).mapIso ((Scheme.Modules.pullback p).mapIso t)
    ≪≫ unitFiberIso p

/-- ★★自明化から得るノルム。 -/
noncomputable def trivNorm (t : L ≅ unitModules V) (p : Spec (CommRingCat.of ℂ) ⟶ V)
    (w : ↥(arcFiber p L)) : ℝ :=
  ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((trivFiberIso t p).hom.hom w)‖

/-- ★★★★切断のノルムは正則関数の絶対値である。 -/
theorem trivNorm_arcEval (t : L ≅ unitModules V) (p : Spec (CommRingCat.of ℂ) ⟶ V)
    (s : (L.val.obj (op ⊤) : Type)) :
    trivNorm t p (arcEval p L s)
      = ‖evalGlobal p ((t.hom.val.app (op ⊤)).hom s)‖ := by
  have h1 : trivNorm t p (arcEval p L s)
      = ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((unitFiberIso p).hom.hom
          (arcEval p (unitModules V) ((t.hom.val.app (op ⊤)).hom s)))‖ :=
    congrArg (fun z => ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((unitFiberIso p).hom.hom z)‖) (arcEval_naturality p t.hom s).symm
  have h2 : ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((unitFiberIso p).hom.hom
        (arcEval p (unitModules V) ((t.hom.val.app (op ⊤)).hom s)))‖
      = ‖evalGlobal p ((t.hom.val.app (op ⊤)).hom s)‖ :=
    congrArg (fun z => ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom z‖)
      (evalUnit_eq p ((t.hom.val.app (op ⊤)).hom s))
  exact h1.trans h2

theorem trivNorm_nonneg (t : L ≅ unitModules V) (p : Spec (CommRingCat.of ℂ) ⟶ V)
    (w : ↥(arcFiber p L)) : 0 ≤ trivNorm t p w :=
  norm_nonneg _

theorem gammaC_hom_inv (c : (CommRingCat.of ℂ : Type)) :
    (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom c) = c :=
  congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) c)
    (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv_hom_id

theorem trivNorm_smul (t : L ≅ unitModules V) (p : Spec (CommRingCat.of ℂ) ⟶ V)
    (c : (CommRingCat.of ℂ : Type)) (w : ↥(arcFiber p L)) :
    trivNorm t p (c • w) = ‖c‖ * trivNorm t p w := by
  have hm : ((trivFiberIso t p).hom).hom (c • w) = c • ((trivFiberIso t p).hom).hom w :=
    map_smul _ c w
  have h0 : (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
        (c • ((trivFiberIso t p).hom).hom w)
      = c * (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (((trivFiberIso t p).hom).hom w) := by
    have h1 := map_mul (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom c)
      (((trivFiberIso t p).hom).hom w)
    rw [gammaC_hom_inv] at h1
    exact h1
  show ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (((trivFiberIso t p).hom).hom (c • w))‖ = _
  rw [hm, h0]
  exact norm_mul _ _

theorem trivNorm_eq_zero_iff (t : L ≅ unitModules V) (p : Spec (CommRingCat.of ℂ) ⟶ V)
    (w : ↥(arcFiber p L)) : trivNorm t p w = 0 ↔ w = 0 := by
  have hb : Function.Bijective (((trivFiberIso t p).hom).hom) :=
    ConcreteCategory.bijective_of_isIso _
  show ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (((trivFiberIso t p).hom).hom w)‖ = 0
    ↔ w = 0
  rw [norm_eq_zero]
  constructor
  · intro h
    refine hb.1 ?_
    have hz : ((trivFiberIso t p).hom).hom (0 : ↥(arcFiber p L)) = 0 := map_zero _
    have hy : ((trivFiberIso t p).hom).hom w = 0 :=
      calc ((trivFiberIso t p).hom).hom w
          = (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom
              ((Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
                (((trivFiberIso t p).hom).hom w)) := (gammaIso_inv_hom _).symm
        _ = (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom 0 := congrArg _ h
        _ = 0 := map_zero _
    rw [hy, hz]
  · intro h
    have hz : ((trivFiberIso t p).hom).hom w = 0 := by
      rw [h]; exact map_zero _
    exact (congrArg (fun y => (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom y) hz).trans
      (map_zero _)

theorem continuous_trivNorm (t : L ≅ unitModules V) (s : (L.val.obj (op ⊤) : Type)) :
    @Continuous _ ℝ (arcTopology V) _ (fun p => trivNorm t p (arcEval p L s)) := by
  letI := arcTopology V
  have he : (fun p : Spec (CommRingCat.of ℂ) ⟶ V => trivNorm t p (arcEval p L s))
      = (fun p => ‖evalGlobal p ((t.hom.val.app (op ⊤)).hom s)‖) :=
    funext fun p => trivNorm_arcEval t p s
  rw [he]
  exact continuous_norm.comp (continuous_evalGlobal _)


/-- ★★★★★★**自明な直線束の計量**——4 法則を束ねたもの。 -/
noncomputable def trivArcMetric (t : L ≅ unitModules V) : ArcMetric V L where
  nrm p := trivNorm t p
  nonneg p := trivNorm_nonneg t p
  eq_zero_iff p := trivNorm_eq_zero_iff t p
  smul p := trivNorm_smul t p

/-! ## ★出典の紐付け(`.src`) -/

def trivArcMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——自明な直線束には連続な計量が入る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
