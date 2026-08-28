/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalRatioCover
import ABC3.Found.GenEll.GlueBump
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★`X_s` の関数は**大域の比**である —— 段 E3c の最後の環（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★これは何か —— 段 E3c の数学が閉じた

`§9-843`（段 E3c 大域版）が要求するのは「試験元 `g` が **`s_j/s_i` の形**である」ことだった。
★本ファイルはそれを**実際に作る**:

> `X_s` の上の関数 `g` に対し、**ある `n`** と `sheafify(M^{⊗(n+1)})` の**大域切断 `t`** があって
>
>     `g = t / (s^{⊗(n+1)})`   （大域の比として）

★★すなわち座標族に `s^{⊗(n+1)}` と `t` を並べれば、`§9-843` の条件が満たされる。

## ★★★★★★組み立て —— 6 本の合流

| 段 | 出典 |
|---|---|
| 分母と重なりに効く単一の指数 | `§9-831` `exists_common_exponent` |
| 指数を 1 つ上げる | `§9-840` `exists_glue_bump_trivValue` |
| 貼った切断の座標は分子 | `§9-839` `trivValue_glued` |
| 被覆で大域の比を同定 | `§9-844` `globalRatio_unique_of_cover` |
| `X_{unit(s^{⊗n})} = X_s` | `§9-844` `nonVanishing_unit_secPow` |
| 比の特徴づけ | `sectionRatio_unique`（在庫、段 D3） |

★**指数を 1 つ上げる**のが要点である——`§9-831` が出す `n` は `0` でありうるが、
`nonVanishing_unit_secPow` は `n ≥ 1` を要る（`n = 0` なら `X_1 = ⊤` になってしまう）。
★★`§9-840` の底上げでちょうど `n+1` になり、しかも分子が `a_i·f_i` になるので
`g·f_i^{n+1} = a_i·f_i` がそのまま出る。

## ★測定の記録

★★★`nonVanishing P (unit s^{⊗n}) ⊓ U_i = X.basicOpen (trivValue M U_i e_i s)` という
**開集合の同一視**が要る（`nonVanishing_inf_unit_secPow`、本ファイル）。
★2 つの開は等しいが**型は違う**ので、`hden` を
`X.presheaf.map (homOfLE (le_of_eq …)).op` で運んでから使う。
★★これが「純粋な配管」と書いていた段の実体である。

★本定理の検査は 34 秒かかる（`maxHeartbeats 4000000`）——
`sheafify(M^{⊗(n+1)})` の項が大きいためである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★開集合の同一視 -/

/-- ★★**`X_{unit(s^{⊗(k+1)})} ⊓ U = D(trivValue s)`**。

★`nonVanishing_inf`（段 D2）＋ `trivValue_unit_secPow`（`§9-839`）＋
`basicOpen_pow_succ`（`§9-844`）の合成である。 -/
theorem nonVanishing_inf_unit_secPow (M : X.PresheafOfModules) (U : X.Opens)
    (e : (restrictPresheafFunctor X U).obj M ≅ 𝟙_ (PresheafModulesOn X U))
    (s : (M.obj (op ⊤) : Type)) (k : ℕ) :
    nonVanishing ((sheafifyFunctor X).obj (presheafTensorPow M (k + 1))).val
        (((sheafifyUnit X (presheafTensorPow M (k + 1))).app (op ⊤)).hom
          (secPow M s (k + 1))) ⊓ U
      = X.basicOpen (trivValue M U e s) := by
  rw [nonVanishing_inf _ U (sheafifyTriv (presheafTensorPow M (k + 1))
    (tensorPowTriv e (k + 1))), trivValue_unit_secPow, basicOpen_pow_succ]

/-! ## ★★★★★★★★★★`X_s` の関数は大域の比である -/

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★**`X_s` の上の関数は、冪を上げれば大域の比になる** —— 段 E3c の最後の環。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`X` が有限個のアフィン開 `U_i`（重なりもアフィン、`M` が自明化）で覆われるとき、
`g ∈ Γ(X, X_s)` に対し ★**ある `n`** と `sheafify(M^{⊗(n+1)})` の**大域切断 `t`** があって

    `g = t / (s^{⊗(n+1)})`   （大域の比 `globalRatio` として）

★★これで `§9-843` の座標条件（試験元が `s_j/s_i` の形）が**実際に満たせる**
——座標族に `s^{⊗(n+1)}` と `t` を並べればよい。

## ★★★機構

`§9-831` で分母と重なりに効く単一の指数 `n` を取り、`§9-840` で `n+1` に底上げする
（★`n = 0` を排除するため。`nonVanishing_unit_secPow` は `n ≥ 1` を要る）。
底上げで分子が `a_i·f_i` になるので `g·f_i^{n+1} = a_i·f_i = trivValue(t)` が出る。
★★あとは `§9-844` の「被覆で一致すれば大域の比」に渡すだけである。 -/
theorem exists_glued_globalRatio {ι : Type} [Fintype ι]
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (g : (Γ(X, nonVanishing M s) : Type)) :
    ∃ (n : ℕ) (t : (((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val.obj
        (op (⊤ : X.Opens)) : Type)),
      X.presheaf.map (homOfLE (le_of_eq (nonVanishing_unit_secPow M hM s n))).op g
        = globalRatio ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
            (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
            t (((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
              (secPow M s (n + 1))) := by
  classical
  obtain ⟨n, a, hden, hagree⟩ := exists_common_exponent Finset.univ M U
    (fun i _ => hU i) (fun i _ j _ => hUij i j) e s
    (fun i => X.presheaf.map (homOfLE (basicOpen_trivValue_le M (U i) (e i) s)).op g)
    (fun i _ j _ => by rw [res_trans, res_trans])
  obtain ⟨t, ht⟩ := exists_glue_bump_trivValue M U hcov e s n 1 a
    (fun i j => hagree i (Finset.mem_univ i) j (Finset.mem_univ j))
  refine ⟨n, t, ?_⟩
  refine globalRatio_unique_of_cover _ _ t _ U
    (fun i => sheafifyTriv (presheafTensorPow M (n + 1)) (tensorPowTriv (e i) (n + 1)))
    (by rw [← inf_iSup_eq, hcov, inf_top_eq]) _ (fun i => ?_)
  refine sectionRatio_unique _ (U i) _ t _ _ ?_
  have hOpen := nonVanishing_inf_unit_secPow M (U i) (e i) s n
  have hd := congrArg (fun z => X.presheaf.map (homOfLE (le_of_eq hOpen)).op z)
    (hden i (Finset.mem_univ i))
  simp only [map_mul, map_pow, algebraMap_basicOpen_eq_res] at hd
  rw [res_trans, res_trans, res_trans] at hd
  rw [trivValue_unit_secPow, ht i, pow_one, map_mul, map_pow, res_trans,
    pow_succ, ← mul_assoc, hd]

/-! ## ★出典の紐付け(`.src`) -/

def nonVanishing_inf_unit_secPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(X_{unit(s^{⊗(k+1)})} ⊓ U = D(trivValue s))",
    sectionId := "genell-prop-1-4" }

def exists_glued_globalRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(X_s の上の関数は冪を上げれば大域の比になる)",
    sectionId := "genell-prop-1-4" }

def exists_glued_globalRatio.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_common_exponent(分母と重なりに効く単一の指数、§9-831)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_common_exponent") 2,
    .citation "[ABC3]" "exists_glue_bump_trivValue(指数を上げる、§9-840)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_glue_bump_trivValue") 2,
    .citation "[ABC3]" "globalRatio_unique_of_cover(被覆で大域の比を同定、§9-844)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_unique_of_cover") 2,
    .citation "[Stacks]" "Lemma 01PW(の消費側——チャートの座標環が生成される段)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)") 7,
    .implicitStep
      ("★**指数を 1 つ上げる**のが要点である——§9-831 が出す n は 0 でありうるが、" ++
       "nonVanishing_unit_secPow は n ≥ 1 を要る(n = 0 なら X_1 = ⊤ になってしまう)。" ++
       "★★§9-840 の底上げでちょうど n+1 になり、しかも分子が a_i·f_i になるので " ++
       "g·f_i^{n+1} = a_i·f_i がそのまま出る") 5,
    .implicitStep
      ("★★★開集合の同一視 nonVanishing P (unit s^{⊗n}) ⊓ U_i = D(trivValue M U_i e_i s) が要る。" ++
       "2 つの開は等しいが**型は違う**ので hden を X.presheaf.map (homOfLE (le_of_eq …)).op で" ++
       "運んでから使う。★これが「純粋な配管」と書いていた段の実体である" ++
       "(検査は 34 秒、maxHeartbeats 4000000)") 4 ]

end ABC3.Found.GenEll
