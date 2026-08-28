/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.SheafifyTrivValue
import ABC3.Meta.Claim

/-!
# ★★★★★★★★層化の側では貼り合わせられる —— 段 E3a-3 の最後の道具（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か —— §9-824 の訂正が要求した段

2026-08-28 に「局所自明な前層加群は層である」が**偽**であると分かった
（反例: `X = {p,q}` 離散、`F({p}) = F({q}) = ℤ`、`F(X) = 0`）。
★★したがって**前層のままでは貼り合わせられない**。大域切断は `sheafify` の側で取る。

★★★本ファイルはその「層化の側で貼る」を型にしたものである:

    重なりで一致する族 `τ_i ∈ (P^sh)(U_i)` は、`⨆ U_i` の上の切断へ**一意に**貼り合う。

## ★★★機構 —— 在庫 2 本の合成

| 道具 | 役割 |
|---|---|
| `((sheafifyFunctor X).obj P).isSheaf` | 層化したものは層である（mathlib） |
| `TopCat.Sheaf.existsUnique_gluing` | 位相空間上の層の貼り合わせ（mathlib） |

★橋渡しは `sheafifyAb`（層化加群層を `TopCat.Sheaf Ab` として読む）だけである。
★★`((sheafifyFunctor X).obj P).val.map` と `(sheafifyAb P).1.map` は**定義的に等しい**
（2026-08-28 実測）ので、加群の言葉のまま `existsUnique_gluing` へ渡せる。

## ★★測定の記録 —— `⊤` へ移す段は `rw` では通らない

`hcov : ⨆ i, U i = ⊤` があっても `rw [← hcov]` は **motive not type correct** で落ちる
——目標に `le_iSup U i : U i ≤ ⨆ U i` という**証明項が現れる**からである（2026-08-28 実測）。
★正しい直し方は「同型で運ぶ」ことである: `⊤ ≤ ⨆ U i` と `⨆ U i ≤ ⊤` の両向きの制限射が
互いに逆であることを使う。`Opens` の射は `Subsingleton` なので合成は `𝟙` に潰れる
（`map_map_apply` / `map_id_apply`、本ファイル）。

## ★残っている段（明示）

★★本ファイルは**貼り合わせの機構だけ**である。段 E3a-3 を閉じるには、これに

1. 有限被覆（段 E2、`exists_finite_cover_of_isAmple`）
2. 分母を払う（段 E3a-1、`exists_pow_mul_eq_res_nonVanishing`）
3. 指数を揃え、重なりで一致させる（§9-826、`exists_common_pow_mul_eq`）
4. 冪と座標の可換性（§9-825、`trivValue_secPow`）

を**組み合わせる**段が要る。★★★4 本とも揃っているので、残りは配管である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★橋渡し —— 層化加群層を `TopCat.Sheaf Ab` として読む -/

/-- ★**層化した加群層の台**を `TopCat.Sheaf Ab X` として読む。

★★これで mathlib の `TopCat.Sheaf.existsUnique_gluing` がそのまま使える。 -/
noncomputable def sheafifyAb (P : X.PresheafOfModules) :
    TopCat.Sheaf Ab X.toPresheafedSpace.carrier :=
  ⟨((sheafifyFunctor X).obj P).presheaf, ((sheafifyFunctor X).obj P).isSheaf⟩

/-! ## ★制限射の関手性（元に当てた形） -/

/-- ★**制限を 2 回するのは 1 回するのと同じ**（元に当てた形）。

★`PresheafOfModules.map` は `restrictScalars` を挟むので `map_comp` が直接使えない
——`presheaf`（`Ab` 値）に降りてから使う。 -/
theorem map_map_apply (M : X.PresheafOfModules) {A B C : X.Opens}
    (f : A ⟶ B) (g : B ⟶ C) (t : (M.obj (op C) : Type)) :
    M.map f.op (M.map g.op t) = M.map (f ≫ g).op t := by
  show M.presheaf.map f.op (M.presheaf.map g.op t) = M.presheaf.map (f ≫ g).op t
  rw [op_comp, M.presheaf.map_comp]
  rfl

/-- ★**同じ開への制限は恒等である**——`Opens` の射は `Subsingleton` だからである。 -/
theorem map_id_apply (M : X.PresheafOfModules) {A : X.Opens} (f : A ⟶ A)
    (t : (M.obj (op A) : Type)) : M.map f.op t = t := by
  show M.presheaf.map f.op t = t
  rw [Subsingleton.elim f (𝟙 A), op_id, M.presheaf.map_id]
  rfl

/-! ## ★★★★★★★★貼り合わせ -/

/-- ★★★★★★★★**層化の側では貼り合わせられる**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

重なり `U i ⊓ U j` で一致する族 `τ i ∈ (P^sh)(U i)` は、
★`⨆ i, U i` の上の切断へ**一意に**貼り合う。

★★これが §9-824 の訂正（「局所自明な前層加群は層である」は偽）が要求した段である
——前層のままでは貼れないので、`sheafify` の側で貼る。 -/
theorem exists_unique_glue {ι : Type} (P : X.PresheafOfModules) (U : ι → X.Opens)
    (τ : ∀ i, (((sheafifyFunctor X).obj P).val.obj (op (U i)) : Type))
    (hcompat : ∀ i j,
      ((sheafifyFunctor X).obj P).val.map
          (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (τ i)
        = ((sheafifyFunctor X).obj P).val.map
          (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (τ j)) :
    ∃! t : (((sheafifyFunctor X).obj P).val.obj (op (⨆ i, U i)) : Type),
      ∀ i, ((sheafifyFunctor X).obj P).val.map (homOfLE (le_iSup U i)).op t = τ i := by
  refine (sheafifyAb P).existsUnique_gluing U τ ?_
  intro i j
  exact hcompat i j

/-- ★★★★★★★★**被覆なら大域切断が出る** —— 段 E3a-3 が求める形。

`⨆ i, U i = ⊤` のとき、重なりで一致する族は
★**`Γ` の側の大域切断**（`(P^sh)(⊤)`）へ一意に貼り合う。

## ★測定の記録

`rw [← hcov]` は **motive not type correct** で落ちる——目標に
`le_iSup U i : U i ≤ ⨆ U i` という証明項が現れるからである（2026-08-28 実測）。
★代わりに両向きの制限射で運ぶ。`Opens` の射は `Subsingleton` なので
合成は `𝟙` に潰れ、`map_id_apply` で消える。 -/
theorem exists_unique_glue_top {ι : Type} (P : X.PresheafOfModules) (U : ι → X.Opens)
    (hcov : (⨆ i, U i) = ⊤)
    (τ : ∀ i, (((sheafifyFunctor X).obj P).val.obj (op (U i)) : Type))
    (hcompat : ∀ i j,
      ((sheafifyFunctor X).obj P).val.map
          (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (τ i)
        = ((sheafifyFunctor X).obj P).val.map
          (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (τ j)) :
    ∃! t : (((sheafifyFunctor X).obj P).val.obj (op (⊤ : X.Opens)) : Type),
      ∀ i, ((sheafifyFunctor X).obj P).val.map (homOfLE (le_top : U i ≤ ⊤)).op t = τ i := by
  obtain ⟨t, ht, hu⟩ := exists_unique_glue P U τ hcompat
  have hup : (⊤ : X.Opens) ≤ ⨆ i, U i := le_of_eq hcov.symm
  refine ⟨((sheafifyFunctor X).obj P).val.map (homOfLE hup).op t, ?_, ?_⟩
  · intro i
    rw [map_map_apply ((sheafifyFunctor X).obj P).val
        (homOfLE (le_top : U i ≤ ⊤)) (homOfLE hup) t,
      Subsingleton.elim (homOfLE (le_top : U i ≤ ⊤) ≫ homOfLE hup) (homOfLE (le_iSup U i))]
    exact ht i
  · intro T' hT'
    have h1 : ((sheafifyFunctor X).obj P).val.map
        (homOfLE (le_top : (⨆ i, U i) ≤ ⊤)).op T' = t := by
      refine hu _ (fun i => ?_)
      rw [map_map_apply ((sheafifyFunctor X).obj P).val
          (homOfLE (le_iSup U i)) (homOfLE (le_top : (⨆ i, U i) ≤ ⊤)) T',
        Subsingleton.elim (homOfLE (le_iSup U i) ≫ homOfLE (le_top : (⨆ i, U i) ≤ ⊤))
          (homOfLE (le_top : U i ≤ ⊤))]
      exact hT' i
    rw [← h1, map_map_apply ((sheafifyFunctor X).obj P).val
      (homOfLE hup) (homOfLE (le_top : (⨆ i, U i) ≤ ⊤)) T',
      map_id_apply ((sheafifyFunctor X).obj P).val _ T']

/-! ## ★出典の紐付け(`.src`) -/

def sheafifyAb.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(層化加群層を TopCat.Sheaf Ab として読む)",
    sectionId := "genell-prop-1-4" }

def exists_unique_glue.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(層化の側では貼り合わせられる)",
    sectionId := "genell-prop-1-4" }

def exists_unique_glue_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(被覆なら大域切断が出る)",
    sectionId := "genell-prop-1-4" }

def map_map_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(制限射の関手性——元に当てた形)",
    sectionId := "genell-prop-1-4" }

def exists_unique_glue.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "TopCat.Sheaf.existsUnique_gluing(位相空間上の層の貼り合わせ)"
      (.inMathlib "TopCat.Sheaf.existsUnique_gluing") 2,
    .citation "[mathlib]" "PresheafOfModules.sheafification(層化したものは層である)"
      (.inMathlib "PresheafOfModules.sheafification") 3,
    .implicitStep
      ("★§9-824 の訂正が要求した段である——「局所自明な前層加群は層である」は**偽**" ++
       "(反例: X = {p,q} 離散、F({p}) = F({q}) = ℤ、F(X) = 0)なので、" ++
       "前層のままでは貼れない。大域切断は sheafify の側で取る") 8,
    .implicitStep
      ("★★rw [← hcov] は motive not type correct で落ちる" ++
       "——目標に le_iSup U i という**証明項が現れる**からである(2026-08-28 実測)。" ++
       "★代わりに両向きの制限射で運ぶ。Opens の射は Subsingleton なので合成は 𝟙 に潰れる") 7,
    .implicitStep
      ("★★★本ファイルは**貼り合わせの機構だけ**である。段 E3a-3 を閉じるには " ++
       "(1) 有限被覆(段 E2)、(2) 分母を払う(§9-822)、(3) 指数を揃え重なりで一致させる(§9-826)、" ++
       "(4) 冪と座標の可換性(§9-825)を組み合わせる段が要る。4 本とも揃っているので残りは配管である") 8 ]

end ABC3.Found.GenEll
