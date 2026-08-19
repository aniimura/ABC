import ABC3.Found.Arakelov.ArcRestrUnit

/-!
# Arakelov (C3) 第 281 ブロック —— **★★★★★★★評価は開集合への制限と両立する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★§9-320 の「残り(1)」が閉じた

    (arcFiberFactor j L p).hom (arcEval (p ≫ j) L s)
      = arcEval p (restrict L j) (restrictSection j L s)

★★★これで **`normSection_continuous` の最後の穴が塞がる**。

## ★★機構——3 段

| 段 | 使うもの | 内容 |
|---|---|---|
| A | 第 278 `unit_conj_compat` + `Adjunction.comp_unit_app` | ★`p ≫ j` の単位を `j` の単位と `p` の単位に分ける |
| B | 第 279 `gamma_arcEval_naturality` | ★`Γ` レベルで層の射と評価を入れ替える |
| C | 第 280 `restrictIso_unit_apply` | ★★`restrictFunctorIsoPullback` を通すと制限写像になる |

★★★`arcFiberFactor` の分解そのもの(`hg`)は **`rfl`** であった
——`moduleSpecΓFunctor` の `map_comp` が `rfl` だからである。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `restrictSection` | ★`⊤` の切断を制限した層の `⊤` の切断へ |
| `unitComp_apply` | ★★段 A(元のレベル) |
| `factor_arcEval` | ★★★★★★★**評価と制限の両立** |

## ★★★測定の記録

§9-320 の見積もり **3–6 ブロック**に対し**実測 4 ブロック**(278–281)。
★★見積もりが当たった理由は §9-324 に書いた通り
——**先に mathlib を探して鍵の在処を確かめてから見積もった**からである。
★★★これは「見積もる前に探す」という手順そのものが精度を上げることを示している。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {U X : Scheme.{0}} (j : U ⟶ X) [IsOpenImmersion j] (L : X.Modules)
  (p : Spec (CommRingCat.of ℂ) ⟶ U)

/-- ★**`⊤` の切断を、制限した層の `⊤` の切断へ運ぶ**。 -/
noncomputable def restrictSection (s : (L.val.obj (op ⊤) : Type)) :
    (((Scheme.Modules.restrict L j).val.obj (op ⊤)) : Type) :=
  (L.val.map (homOfLE (le_top : (j ''ᵁ (⊤ : U.Opens)) ≤ ⊤)).op).hom s

/-- ★★**段 A**——`p ≫ j` での評価を、`pullbackComp` を通して
「`j` の単位のあと `p` で評価」に分ける。 -/
theorem unitComp_apply (s : (L.val.obj (op ⊤) : Type)) :
    (((Scheme.Modules.pullbackComp p j).inv.app L).val.app (op ⊤)).hom (arcEval (p ≫ j) L s)
      = arcEval p ((Scheme.Modules.pullback j).obj L)
        ((((Scheme.Modules.pullbackPushforwardAdjunction j).unit.app L).val.app (op ⊤)).hom s) := by
  have h := congrArg
    (fun (m : L ⟶ (Scheme.Modules.pushforward (p ≫ j)).obj
        ((Scheme.Modules.pullback j ⋙ Scheme.Modules.pullback p).obj L)) =>
      (m.val.app (op ⊤)).hom s) (unit_conj_compat p j L)
  simp only [Adjunction.comp_unit_app] at h
  exact h.symm

/-- ★★★★★★★**評価は開集合への制限と両立する**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★「`X` の切断を点 `p ≫ j` で評価する」ことと
「`U` へ制限してから点 `p` で評価する」ことが、
ファイバーの分解同型 `arcFiberFactor` の下で一致する。

★★これが計量の貼り合わせで **`X` 側の切断と `U` 側のノルムを繋ぐ**唯一の橋である。 -/
theorem factor_arcEval (s : (L.val.obj (op ⊤) : Type)) :
    (arcFiberFactor j L p).hom.hom (arcEval (p ≫ j) L s)
      = arcEval p (Scheme.Modules.restrict L j) (restrictSection j L s) := by
  have hg : (arcFiberFactor j L p).hom.hom (arcEval (p ≫ j) L s)
      = ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).map
          ((Scheme.Modules.pullback p).map
            (((Scheme.Modules.restrictFunctorIsoPullback j).app L).inv))).hom
        ((((Scheme.Modules.pullbackComp p j).inv.app L).val.app (op ⊤)).hom
          (arcEval (p ≫ j) L s)) := rfl
  rw [hg, unitComp_apply, gamma_arcEval_naturality, restrictIso_unit_apply]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def factor_arcEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——評価が開集合への制限と両立すること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
