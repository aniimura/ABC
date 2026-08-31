/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ChartSurjectiveCoords
import ABC3.Found.Arakelov.SheafifyTrivValue
import ABC3.Found.GenEll.Globalize
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★貼った切断の座標は分子そのもの —— 段 E3c の座標条件を作る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か —— 環が閉じた

`§9-838` は段 E3c を「**座標条件** `g·s_i = s_j` が満たせるか」に落とした。
★本ファイルは**その条件を実際に作る**:

> `X_s` の上の関数 `g` に対し、**ある `n`** と `sheafify(M^{⊗n})` の**大域切断 `t`** があって、
> どのチャートでも `g · trivValue(s^{⊗n}) = trivValue(t)`。

★★すなわち `t` を座標の 1 つに取れば、`§9-837` により `g` は `A⁰_{x_i}` の像に入る。

## ★★★★★★機構 —— 3 本の合成

| 段 | 出典 | 役割 |
|---|---|---|
| 大域化 | `§9-831` `exists_global_section_of_localData` | `t` と分子 `a_i` を作る |
| ★**貼った切断の座標は分子** | 本ファイル `trivValue_glued` | `trivValue(t) = a_i` |
| ★★冪と層化 | `§9-825` ＋ `§9-823` | `trivValue(unit s^{⊗n}) = (trivValue s)^n` |

★`trivValue_glued` が本ファイルの新しい部分である。
その中身は **`trivEquiv_sheafifyTriv`（切断版）** ——`§9-823` は大域切断についてだったが、
貼り合わせで出てくるのは**チャート上の切断**なので、切断版が要る。
★★証明は `§9-823` と同型（`trivEquiv_isoComp` ＋ `hom_inv_id`）である。

## ★残っている段（明示）

★★★**残るのは有限の帳簿だけ**である:

1. 試験元は有限個（`§9-832`）、チャートも有限個（`§9-817`）なので、
   得られる `t` も有限個——それらを 1 つの座標族 `s' : Fin (N+1) → …` にまとめる
2. ★ただし `g` ごとに `n` が違う。**指数を揃える段が要る**
   ——`t ↦ t ⊗ (s^{⊗k})` で `n` を上げる操作を作ればよい（`§9-826` の型の議論）

★★1 は純粋な帳簿だが、2 は**まだ書いていない**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

/-! ## ★★★★★層化の単位は座標を変えない（切断版） -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★**層化の単位は座標を変えない**（切断版）。

`§9-823` の `trivValue_sheafifyTriv` は**大域切断**についてだったが、
★貼り合わせで出てくるのは**チャート上の切断**なので、こちらの形が要る。
★★証明は同型である——`trivEquiv_isoComp` ＋ `hom_inv_id`。 -/
theorem trivEquiv_sheafifyTriv {X : Scheme.{0}} (P : X.PresheafOfModules) {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V))
    (x : (P.obj (op V) : Type)) :
    trivEquiv ((sheafifyFunctor X).obj P).val V (sheafifyTriv P e)
        (((sheafifyUnit X P).app (op V)).hom x)
      = trivEquiv P V e x := by
  haveI inst := isIso_restrictMap_sheafifyUnit X V P (isSheaf_restrict_of_triv X V P e)
  rw [sheafifyTriv_eq_of P inst e]
  show trivEquiv _ V ((@asIso _ _ _ _ _ inst).symm ≪≫ e) _ = _
  rw [trivEquiv_isoComp]
  congr 1
  exact congrArg (fun (χ : (restrictPresheafFunctor X V).obj P
      ⟶ (restrictPresheafFunctor X V).obj P) => χ.app (op (Over.mk (𝟙 V))) x)
    (@asIso _ _ _ _ _ inst).hom_inv_id

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★貼った切断の座標は分子そのもの -/

/-- ★★★★★★★★**貼った大域切断のチャート上の座標は、分母を払って得た分子そのものである**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★これで `§9-831`（大域化）の出力が `§9-837`（比の判定）の入力に**そのまま入る**。 -/
theorem trivValue_glued (M : X.PresheafOfModules) (U : X.Opens) (n : ℕ)
    (e : (restrictPresheafFunctor X U).obj M ≅ 𝟙_ (PresheafModulesOn X U))
    (a : (Γ(X, U) : Type))
    (t : (((sheafifyFunctor X).obj (presheafTensorPow M n)).val.obj
      (op (⊤ : X.Opens)) : Type))
    (ht : ((sheafifyFunctor X).obj (presheafTensorPow M n)).val.map
        (homOfLE (le_top : U ≤ ⊤)).op t
      = ((sheafifyUnit X (presheafTensorPow M n)).app (op U)).hom
          (secOfFun M U e n a)) :
    trivValue ((sheafifyFunctor X).obj (presheafTensorPow M n)).val U
        (sheafifyTriv (presheafTensorPow M n) (tensorPowTriv e n)) t = a := by
  rw [trivValue_eq_trivEquiv]
  show trivEquiv _ U _ (((sheafifyFunctor X).obj (presheafTensorPow M n)).val.map
      (homOfLE (le_top : U ≤ ⊤)).op t) = a
  rw [ht, trivEquiv_sheafifyTriv, trivEquiv_secOfFun]

/-- ★★★★**分母の側の座標は `(trivValue s)^n` である**（`§9-823` ＋ `§9-825`）。 -/
theorem trivValue_unit_secPow (M : X.PresheafOfModules) (U : X.Opens) (n : ℕ)
    (e : (restrictPresheafFunctor X U).obj M ≅ 𝟙_ (PresheafModulesOn X U))
    (s : (M.obj (op ⊤) : Type)) :
    trivValue ((sheafifyFunctor X).obj (presheafTensorPow M n)).val U
        (sheafifyTriv (presheafTensorPow M n) (tensorPowTriv e n))
        (((sheafifyUnit X (presheafTensorPow M n)).app (op ⊤)).hom (secPow M s n))
      = (trivValue M U e s) ^ n := by
  rw [trivValue_sheafifyTriv, trivValue_secPow]

/-! ## ★★★★★★★★★★段 E3c の座標条件を実際に作る -/

/-- ★★★★★★★★★★**段 E3c の座標条件を実際に作る** —— 環が閉じた。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`X_s` の上の関数 `g` に対し、★**ある `n`** と `sheafify(M^{⊗n})` の**大域切断 `t`** があって、
どのチャート `U_i` でも

    `g|·(trivValue(unit s^{⊗n}))| = (trivValue t)|`   （`D(f_i)` の上で）

★★すなわち `t` を座標の 1 つに取れば、`§9-837` により `g` は `A⁰_{x_i}` の像に入る
——`§9-838` の座標条件が満たされる。

★★★これで **`§9-822`→`§9-826`→`§9-827`→`§9-828`→`§9-829`→`§9-831`→`§9-837`→`§9-838`**
の環が閉じた。残るのは**有限の帳簿と指数揃え**だけである。 -/
theorem exists_glued_ratio_section {ι : Type} [Fintype ι]
    (M : X.PresheafOfModules) (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (g : (Γ(X, nonVanishing M s) : Type)) :
    ∃ (n : ℕ) (t : (((sheafifyFunctor X).obj (presheafTensorPow M n)).val.obj
        (op (⊤ : X.Opens)) : Type)), ∀ i,
      X.presheaf.map (homOfLE (basicOpen_trivValue_le M (U i) (e i) s)).op g
          * X.presheaf.map (homOfLE (X.basicOpen_le (trivValue M (U i) (e i) s))).op
            (trivValue ((sheafifyFunctor X).obj (presheafTensorPow M n)).val (U i)
              (sheafifyTriv (presheafTensorPow M n) (tensorPowTriv (e i) n))
              (((sheafifyUnit X (presheafTensorPow M n)).app (op ⊤)).hom (secPow M s n)))
        = X.presheaf.map (homOfLE (X.basicOpen_le (trivValue M (U i) (e i) s))).op
            (trivValue ((sheafifyFunctor X).obj (presheafTensorPow M n)).val (U i)
              (sheafifyTriv (presheafTensorPow M n) (tensorPowTriv (e i) n)) t) := by
  obtain ⟨n, a, t, hden, ht⟩ :=
    exists_global_section_of_globalFun M U hcov hU hUij e s g
  refine ⟨n, t, fun i => ?_⟩
  rw [trivValue_glued M (U i) n (e i) (a i) t (ht i), trivValue_unit_secPow, map_pow]
  exact hden i

/-! ## ★出典の紐付け(`.src`) -/

def trivEquiv_sheafifyTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(層化の単位は座標を変えない——切断版)",
    sectionId := "genell-prop-1-4" }

def trivValue_glued.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(貼った大域切断の座標は分子そのもの)",
    sectionId := "genell-prop-1-4" }

def trivValue_unit_secPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(分母の側の座標は (trivValue s)^n)",
    sectionId := "genell-prop-1-4" }

def exists_glued_ratio_section.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3c の座標条件を実際に作る)",
    sectionId := "genell-prop-1-4" }

def exists_glued_ratio_section.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_global_section_of_localData(段 E3a-3、§9-831)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_global_section_of_localData") 2,
    .citation "[ABC3]" "trivValue_secPow(冪と座標は可換、§9-825)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_secPow") 2,
    .citation "[ABC3]" "trivValue_sheafifyTriv(座標は層化で変わらない、§9-823)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivValue_sheafifyTriv") 2,
    .citation "[Stacks]" "Lemma 01PW(の消費側)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)") 7,
    .implicitStep
      ("★★★**残るのは有限の帳簿と指数揃えだけ**である: " ++
       "(1) 試験元は有限個(§9-832)、チャートも有限個(§9-817)なので得られる t も有限個" ++
       "——それらを 1 つの座標族にまとめる(純粋な帳簿)。" ++
       "(2) ★ただし g ごとに n が違う。**指数を揃える段が要る**" ++
       "——t ↦ t ⊗ (s^{⊗k}) で n を上げる操作を作ればよい(§9-826 の型の議論)。" ++
       "★★(2) はまだ書いていない") 7 ]

end ABC3.Found.GenEll
