/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.CommonDegree
import ABC3.Found.Arakelov.TrivTensor
import ABC3.Meta.Claim

/-!
# ★★★★★テンソル冪の自明化と座標 `trivValue (s^{⊗n}) = (trivValue s)^n`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★これは段 E3a-3 の道具である

段 E3a-3（大域化）では、各チャート `U_j` で得た**関数** `a_j ∈ Γ(X, U_j)` を
`M^{⊗n}` の切断に戻して貼り合わせる。★その戻し方が

    `(trivEquiv (M^{⊗n}) U_j (tensorPowTriv e_j n)).symm a_j`

であり、`tensorPowTriv`（本ファイル）がその自明化を与える。
★★そして**重なりで一致するか**を見るには、両側の座標を比べる必要がある
——そこで効くのが `trivValue_secPow`（本ファイル）である:

    `trivValue (M^{⊗n}) V (tensorPowTriv e n) (s^{⊗n}) = (trivValue M V e s)^n`

★★★すなわち**`n` 冪を取ることと座標を読むことは可換である**。
遷移単元（段 D1、`trivValue_congr`）が `n` 乗で効く、というのが貼り合わせの機構になる。

## ★★★機構 —— 2 本だけ

| 段 | 道具 |
|---|---|
| 後継 `n+1` | `trivValue_tensor`（`TrivTensor.lean`、`rfl`）——座標は掛け算に落ちる |
| 基点 `n = 0` | `trivValue_unit'`（本ファイル）——単位前層の座標は**ただの制限**である |

`presheafTensorPow M (n+1) = M ⊗ presheafTensorPow M n` は `rfl` なので、
`tensorTriv` を再帰で重ねるだけで自明化ができる。

## ★測定の記録（`tools/lean-idioms.md` にも入れた）

★★**`restrictPresheafUnit` の向きは `𝟙_ ≅ (restrict).obj (𝟙_)` である**
——「`(restrict).obj (𝟙_) ≅ 𝟙_`」ではない。★基点を `restrictPresheafUnit` と書いて
`trivValue … = 1` が `rfl` でも `simp` でも落ちず半日詰まったが、原因は向きだけであった
（2026-08-28 実測）。★★★在庫の `trivValue_unit_top` が `.symm` を使っていたのが合図であり、
**そこを読んでいれば即座に分かった**——CLAUDE.md の在庫の規律どおり、まず引くこと。

★また `trivValue_unit_top`（在庫、`AmpleDef.lean`）は `V = ⊤` に限られていた。
本ファイルの `trivValue_unit'` はそれを**任意の `V`** へ広げたものである
（`V = ⊤` なら制限が恒等なので在庫と一致する）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★テンソル冪の自明化 -/

/-- ★★★★**テンソル冪の自明化** `(M^{⊗n})|_V ≅ 𝟙_`。

★`presheafTensorPow M (n+1) = M ⊗ presheafTensorPow M n` は `rfl` なので、
`tensorTriv` を再帰で重ねるだけである。

★★基点は `restrictPresheafUnit.symm` である——`restrictPresheafUnit` の向きは
`𝟙_ ≅ (restrict).obj (𝟙_)` なので、**`.symm` が要る**。 -/
noncomputable def tensorPowTriv {M : X.PresheafOfModules} {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    ∀ n : ℕ, (restrictPresheafFunctor X V).obj (presheafTensorPow M n)
      ≅ 𝟙_ (PresheafModulesOn X V)
  | 0 => restrictPresheafUnit.symm
  | n + 1 => tensorTriv e (tensorPowTriv e n)

theorem tensorPowTriv_zero {M : X.PresheafOfModules} {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    tensorPowTriv e 0 = restrictPresheafUnit.symm := rfl

theorem tensorPowTriv_succ {M : X.PresheafOfModules} {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) (n : ℕ) :
    tensorPowTriv e (n + 1) = tensorTriv e (tensorPowTriv e n) := rfl

/-! ## ★★★★★単位前層の座標は「ただの制限」である -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★**単位前層の座標はただの制限である**——任意の開 `V` について。

    `trivValue (𝟙_) V restrictPresheafUnit.symm s = s|_V`

★在庫の `trivValue_unit_top`（`AmpleDef.lean`）は `V = ⊤` に限られていた。
本補題はそれを任意の `V` へ広げたものである（`V = ⊤` なら制限が恒等）。 -/
theorem trivValue_unit' (V : X.Opens)
    (s : (𝟙_ X.PresheafOfModules).obj (op (⊤ : X.Opens))) :
    trivValue (𝟙_ X.PresheafOfModules) V restrictPresheafUnit.symm s
      = X.presheaf.map (homOfLE (le_top : V ≤ ⊤)).op s := by
  simp [trivValue, secOn, restrictPresheafUnit, Functor.Monoidal.εIso]
  rfl

/-- ★★★★**したがって単位切断 `1` の座標は `1` である**——制限は環準同型だからである。 -/
theorem trivValue_unit_one (V : X.Opens) :
    trivValue (𝟙_ X.PresheafOfModules) V restrictPresheafUnit.symm
      ((1 : (Γ(X, (⊤ : X.Opens)) : Type))) = 1 := by
  rw [trivValue_unit' V]
  exact map_one (X.presheaf.map (homOfLE (le_top : V ≤ ⊤)).op).hom

/-! ## ★★★★★★★★冪を取ることと座標を読むことは可換である -/

/-- ★★★★★**後継段——座標は掛け算に落ちる**。

    `trivValue (M^{⊗(n+1)}) (s^{⊗(n+1)})
       = trivValue M s · trivValue (M^{⊗n}) (s^{⊗n})`

★機構は `trivValue_tensor`（`TrivTensor.lean`、`rfl`）だけである。 -/
theorem trivValue_secPow_succ {M : X.PresheafOfModules} {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : (M.obj (op ⊤) : Type)) (n : ℕ) :
    trivValue (presheafTensorPow M (n + 1)) V (tensorPowTriv e (n + 1)) (secPow M s (n + 1))
      = trivValue M V e s
        * trivValue (presheafTensorPow M n) V (tensorPowTriv e n) (secPow M s n) := by
  show trivValue (M ⊗ presheafTensorPow M n) V (tensorTriv e (tensorPowTriv e n))
    (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] secPow M s n) = _
  rw [trivValue_tensor]

/-- ★★★★★★★★**冪を取ることと座標を読むことは可換である** —— 段 E3a-3 の鍵。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `trivValue (M^{⊗n}) V (tensorPowTriv e n) (s^{⊗n}) = (trivValue M V e s)^n`

★★これで「`n` 冪へ移ってから貼り合わせる」議論が**座標の側で `n` 乗するだけ**になる。
遷移単元（段 D1、`trivValue_congr`）も `n` 乗で効くので、重なりの一致が見える。

★★★機構は 2 本だけ——後継は `trivValue_tensor`（`rfl`）、基点は `trivValue_unit_one`。 -/
theorem trivValue_secPow {M : X.PresheafOfModules} {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : (M.obj (op ⊤) : Type)) : ∀ n : ℕ,
      trivValue (presheafTensorPow M n) V (tensorPowTriv e n) (secPow M s n)
        = (trivValue M V e s) ^ n
  | 0 => by
      show trivValue (𝟙_ X.PresheafOfModules) V restrictPresheafUnit.symm
        ((1 : (Γ(X, (⊤ : X.Opens)) : Type))) = _
      rw [pow_zero, trivValue_unit_one]
  | n + 1 => by
      show trivValue (M ⊗ presheafTensorPow M n) V (tensorTriv e (tensorPowTriv e n))
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] secPow M s n) = _
      rw [trivValue_tensor, trivValue_secPow e s n, pow_succ]
      ring

/-! ## ★出典の紐付け(`.src`) -/

def tensorPowTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(語彙——テンソル冪の自明化)",
    sectionId := "genell-prop-1-4" }

def trivValue_unit'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(単位前層の座標はただの制限である)",
    sectionId := "genell-prop-1-4" }

def trivValue_secPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(冪を取ることと座標を読むことは可換)",
    sectionId := "genell-prop-1-4" }

def trivValue_secPow.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivValue_tensor(座標の読みは掛け算に落ちる、TrivTensor.lean)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivValue_tensor") 6,
    .citation "[ABC3]" "trivValue_unit_top(V = ⊤ の場合、AmpleDef.lean)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_unit_top") 6,
    .implicitStep
      ("★restrictPresheafUnit の向きは 𝟙_ ≅ (restrict).obj (𝟙_) である" ++
       "——逆ではない。基点に .symm が要る(2026-08-28 実測)。" ++
       "★★在庫の trivValue_unit_top が .symm を使っていたのが合図であった") 7,
    .implicitStep
      ("★★段 E3a-3(大域化)で使う: 各チャートの関数 a_j を " ++
       "(trivEquiv (M^{⊗n}) U_j (tensorPowTriv e_j n)).symm で切断に戻し、" ++
       "重なりの一致を座標の n 乗で見る。★貼り合わせは §9-824 により " ++
       "sheafify (M^{⊗n}) の側で行う") 8 ]

end ABC3.Found.GenEll
