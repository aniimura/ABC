/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.SheafifyTriv
import ABC3.Found.GenEll.AmpleDef
import ABC3.Meta.Claim

/-!
# ★★★★★★`trivValue` は層化で変わらない（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か —— 段 E3a-2 の本体

`ResearchPaper/mathlib-gap.json` の `ample-and-projective-embedding` の測定（2026-08-28）:

* 「局所自明な前層加群は層である」は**偽**である（反例あり）。
* したがって大域切断は**層化 `sheafify M` の側**で取る。
* そのとき**局所の理論が層化で変わらない**ことが要る。

★本ファイルはその核である:

    `trivValue (P^sh) V (sheafifyTriv P e) (unit s) = trivValue P V e s`

★★これがあれば `nonVanishing`・`homogValue`・`homogRatio` はすべて層化の側へ移る
——どれも `trivValue` で書かれているからである。

## ★★★★★機構 —— 在庫 3 本の合成

| 道具 | 役割 |
|---|---|
| `sheafifyTriv`（`SheafifyTriv.lean`） | 層化した側の自明化 `= (制限した単位の asIso).symm ≪≫ e` |
| `trivEquiv_isoComp`（同上） | 自明化に同型を前置すると座標が前置分ずれる（`rfl`） |
| `PresheafOfModules.naturality_apply`（mathlib） | 単位の自然性——`secOn` が単位と可換 |

★あとは `hom_inv_id` で `asIso` を潰すだけである。

## ★測定の記録

`PresheafOfModules.naturality_apply` の向きは
`f.app V (P.map g x) = Q.map g (f.app ⊤ x)` である。★`secOn` の形に合わせるには
**`.symm` を取る**（2026-08-28 実測）。

## ★残っている段（明示）

★★`nonVanishing (P^sh) (unit s) = nonVanishing P s` は本ファイルに無い。
`§9-818` の `nonVanishing_tensor` と**同じ型**の議論で出る
——「各点で `P` が自明化する開を取り、そこで両辺を `basicOpen` に直す」。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.GenEll

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**`trivValue` は層化で変わらない**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `trivValue (P^sh) V (sheafifyTriv P e) (unit s) = trivValue P V e s`

★★これで `nonVanishing`・`homogValue`・`homogRatio` が層化の側へ移せる
——どれも `trivValue` で書かれているからである。

★★★機構は `trivEquiv_isoComp`（前置した同型の分だけ座標がずれる）＋
単位の自然性（`secOn` が単位と可換）＋ `hom_inv_id` だけである。 -/
theorem trivValue_sheafifyTriv {X : Scheme.{0}} (P : X.PresheafOfModules) {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V))
    (s : (P.obj (op ⊤) : Type)) :
    trivValue ((sheafifyFunctor X).obj P).val V (sheafifyTriv P e)
        (((sheafifyUnit X P).app (op ⊤)).hom s)
      = trivValue P V e s := by
  haveI inst := isIso_restrictMap_sheafifyUnit X V P (isSheaf_restrict_of_triv X V P e)
  rw [sheafifyTriv_eq_of P inst e, trivValue_eq_trivEquiv, trivValue_eq_trivEquiv]
  show trivEquiv _ V ((@asIso _ _ _ _ _ inst).symm ≪≫ e) _ = _
  rw [trivEquiv_isoComp]
  congr 1
  have hnat := PresheafOfModules.naturality_apply (sheafifyUnit X P)
    (homOfLE (le_top : V ≤ ⊤)).op s
  have hs : secOn ((sheafifyFunctor X).obj P).val V (((sheafifyUnit X P).app (op ⊤)).hom s)
      = ((restrictPresheafFunctor X V).map (sheafifyUnit X P)).app (op (Over.mk (𝟙 V)))
          (secOn P V s) := hnat.symm
  rw [hs]
  exact congrArg (fun (χ : (restrictPresheafFunctor X V).obj P
      ⟶ (restrictPresheafFunctor X V).obj P) => χ.app (op (Over.mk (𝟙 V))) (secOn P V s))
    (@asIso _ _ _ _ _ inst).hom_inv_id

/-- ★★★**したがって `basicOpen` も変わらない**。 -/
theorem basicOpen_trivValue_sheafifyTriv {X : Scheme.{0}} (P : X.PresheafOfModules)
    {V : X.Opens} (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V))
    (s : (P.obj (op ⊤) : Type)) :
    X.basicOpen (trivValue ((sheafifyFunctor X).obj P).val V (sheafifyTriv P e)
        (((sheafifyUnit X P).app (op ⊤)).hom s))
      = X.basicOpen (trivValue P V e s) := by
  rw [trivValue_sheafifyTriv]

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**非消失軌跡は層化で変わらない** —— 段 E3a-2 が閉じた。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `X_{unit s}` （層化の側） `= X_s` （前層の側）

★★これで `§9-817`〜`§9-820`（被覆・次数揃え）と `§9-822`（分母を払う）が
**層化の側へそのまま移る**。

★★★機構は `§9-818` の `nonVanishing_tensor` と同じ型である
——「各点で `P` が自明化する開 `W` を取り、`V ⊓ W` の上で両辺を `basicOpen` に直す」。 -/
theorem nonVanishing_sheafify {X : Scheme.{0}} (P : X.PresheafOfModules)
    (hP : IsLocallyTrivial X P) (s : (P.obj (op ⊤) : Type)) :
    nonVanishing ((sheafifyFunctor X).obj P).val
        (((sheafifyUnit X P).app (op ⊤)).hom s)
      = nonVanishing P s := by
  apply le_antisymm
  · intro x hx
    obtain ⟨V, e', hxe⟩ := (mem_nonVanishing_iff _ _ x).1 hx
    obtain ⟨S, hS, htriv⟩ := hP ⊤
    obtain ⟨W, i, hi, hxW⟩ := hS x trivial
    obtain ⟨e⟩ := htriv i hi
    have hxV : x ∈ V := X.basicOpen_le _ hxe
    have hmem : x ∈ nonVanishing ((sheafifyFunctor X).obj P).val
        (((sheafifyUnit X P).app (op ⊤)).hom s) ⊓ (V ⊓ W) := ⟨hx, hxV, hxW⟩
    rw [nonVanishing_inf _ (V ⊓ W)
      (sheafifyTriv P (trivialOfLe P (inf_le_right : V ⊓ W ≤ W) e)) _,
      trivValue_sheafifyTriv,
      ← nonVanishing_inf P (V ⊓ W) (trivialOfLe P (inf_le_right : V ⊓ W ≤ W) e) s] at hmem
    exact hmem.1
  · intro x hx
    obtain ⟨W, e, hxe⟩ := (mem_nonVanishing_iff P s x).1 hx
    refine basicOpen_trivValue_le ((sheafifyFunctor X).obj P).val W (sheafifyTriv P e) _ ?_
    rw [trivValue_sheafifyTriv]
    exact hxe

/-! ## ★出典の紐付け(`.src`) -/

def nonVanishing_sheafify.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(非消失軌跡は層化で変わらない)",
    sectionId := "genell-prop-1-4" }

def trivValue_sheafifyTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(trivValue は層化で変わらない)",
    sectionId := "genell-prop-1-4" }

def basicOpen_trivValue_sheafifyTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(したがって basicOpen も層化で変わらない)",
    sectionId := "genell-prop-1-4" }

def trivValue_sheafifyTriv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "sheafifyTriv / trivEquiv_isoComp(SheafifyTriv.lean)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.sheafifyTriv") 3,
    .citation "[mathlib]" "PresheafOfModules.naturality_apply(単位の自然性)"
      (.inMathlib "PresheafOfModules.naturality_apply") 7,
    .implicitStep
      ("★PresheafOfModules.naturality_apply の向きは f.app V (P.map g x) = Q.map g (f.app ⊤ x) " ++
       "である。secOn の形に合わせるには **.symm を取る**(2026-08-28 実測)") 7,
    .implicitStep
      ("★★nonVanishing (P^sh) (unit s) = nonVanishing P s は本ファイルに無い。" ++
       "§9-818 の nonVanishing_tensor と**同じ型**の議論で出る" ++
       "——「各点で P が自明化する開を取り、そこで両辺を basicOpen に直す」") 7 ]

end ABC3.Found.Arakelov
