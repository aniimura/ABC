import ABC3.Found.PGC.GainedToWildStep
import ABC3.Found.PGC.PureStepSetup

/-!
# [pGC] ★★「26 個の仮説」の数え方が誤り —— 正しくは 13、そして深さ 1 だけでは足りない

## 持ち場（前波で「次の 1 点」とした点）

前波の私の言葉（逐語）:

> `Gained*` 族のどれか 1 件を `F := K.carrier`、`M := K(x)`（`[M:K] = p`）に代入し、
> ★その **26 個の仮説**のうち何個が `PAdicLocalField` の設定で自動的に満たされるかを数えること。

## ★★★訂正 —— 問いの立て方が誤りだった（自己訂正 18 度目）

★**数えるべきは 26 ではなく 13。**
`PureStepSetup.lean:309` の `structure WildStep` が**すでに畳んでいる**。フィールドは

  データ 4: `unif`, `gen`, `i`, `e`
  証明 9 : `hi`, `he`, `hdeg`, `htop`, `hiso`, `hjump`, `hnormp`, `heM`, `hval`

の**計 13**（`sed -n '309,328p' PureStepSetup.lean` で数えた）。
★同ファイルの docstring も「`GainedJumpFree` の配管 10 本を 1 つの構造にまとめ、
うち `hsep`/`hsp`/`hconj`/`hne` は**定理として出す**」と言っている。
⇒ ★私が数えた「26」は `GainedDescentBridge.lean:527` の**畳む前**の姿だった。

★★どちらを取るかの判断（持ち場で私が決めることになっていた点）:
★**`PureStepSetup` を取る**。理由は 2 つ、いずれも測った:

1. 結論が `∏` ではなく ★**`axDecay p 1` そのもの**（`:414`
   `exists_norm_sub_algebraMap_le_axDecay : ‖x − y‖ ≤ axDecay p 1 * ‖S.gen x − x‖`）で、
   ★私の `GainedToWildStep.axWildDescent_of_depth_one_approx` の `C` にそのまま入る。
2. 仮説が 13 に畳まれており、うち ★**2 つは同ファイル §5 が `PAdicLocalField` から供給済み**:
   `hnormp` ← `norm_natCast_p_eq_inv`（`:455`）、
   `hiso` ← `norm_algHom_eq_of_pAdicLocalField`（`:466`）。

## ★測定（コマンドと結果）

| 測ったもの | 結果 |
|---|---|
| `WildStep` のフィールド数 | ★**13**（データ 4 ＋ 証明 9） |
| `WildStep` を**構成する**側 | ★★**0 件**（`grep -rn "WildStep" lean/ABC3/Found/PGC/*.lean \| grep -cv PureStepSetup.lean` → `2`、いずれも本ファイル群の `namespace` 行） |
| `PAdicLocalField` から供給済みのフィールド | ★**2 / 13**（`hnormp`, `hiso`） |
| `hval` の形（`∀ c ≠ 0, ∃ m : ℤ, ‖algebraMap c‖ = ‖π‖^(p·m)`）を持つファイル | **39**（`grep -rln "∃ m : ℤ, ‖algebraMap" ... \| wc -l`）。代表 `TotallyRamifiedValueGroup.lean:276 exists_zpow_norm_intermediate` |

★★**`WildStep` は木の中で一度も構成されていない。** 仕様は在るが実例が無い。

## ★★★本波の決定的な測定 —— 深さ 1 だけでは足りない（§2、証明つき）

`PureStepSetup` の道は ★**`k = 0` 専用**である（同ファイル `:430` の docstring:
「`k ≥ 1` では `hne` が偽なので、これが `_of_conj_pure` を非空虚に使える**唯一の**場合」）。
⇒ 繋げられるのは **wild 深さ 1 の段だけ**で、定数は `axDecay p 1`、他の深さは `p` のまま。

★★そこで測った: `c k = if k = 1 then axDecay p 1 else p` の有限積は ★**非有界**（§2）。
⇒ ★★**深さ 1 を `p` から `axDecay p 1` に落としても `AxLemma` は出ない。**
★**すべての深さ `k` で `axDecay p k` が要る**ことが、これで型になった。

★これは `WildDepthFieldDescent.lean` の docstring
（`p^n` は非有界なので `AxLemma` は出ない）と同じ形の主張だが、
★**「深さ 1 だけ改善した場合」については木にも私にも無かった**。

## 逸脱の記録

- §1 は `GainedToWildStep.axWildDescent_of_depth_one_approx` に `C := axDecay p 1` を
  入れるだけで、★新しい数学は 0。
- §2 は木の `AxTowerDecay.prod_unbounded_of_one_lt`（`:272`）に
  `D := min (axDecay p 1) p` を入れた。★`D > 1` は `axDecay p 1 = p^{1/(p−1)} > 1` から。
- ★本ファイルは `WildStep` の**実例を作っていない**（13 フィールドのうち 11 が未供給）。
  ★作れると書いていない。★測っただけである。
-/

namespace ABC3.Found.PGC

namespace WildStepFieldSupply

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 深さ 1 の定数を `axDecay p 1` にする（`PureStepSetup` の結論の形） -/

section DepthOne

/-- ★`PureStepSetup.WildStep.exists_norm_sub_algebraMap_le_axDecay`（`:414`）の結論を
`K.closure` の深さ 1 の元について仮定すると、深さ 1 の定数が `axDecay p 1` になる。

★他の深さは在庫の `axWildDescent_prime`（定数 `p`、無条件）のまま。 -/
theorem axWildDescent_of_wildStep_approx (K : PAdicLocalField p)
    (happ : ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      wildDepth K x = 1 → ∃ a : K.carrier,
        ‖x - algebraMap K.carrier K.closure a‖ ≤ axDecay p 1 * ε) :
    AxWildDescent K (fun k => if k = 1 then axDecay p 1 else (p : ℝ)) :=
  GainedToWildStep.axWildDescent_of_depth_one_approx K
    (le_trans zero_le_one (one_le_axDecay p 1)) happ

end DepthOne

/-! ## §2 ★★★深さ 1 だけでは `AxLemma` は出ない（非有界性の証明） -/

section NotEnough

/-- `1 < axDecay p 1`（`= p^{1/(p−1)}`）。 -/
theorem one_lt_axDecay_one : (1 : ℝ) < axDecay p 1 := by
  have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  rw [axDecay_one]
  exact Real.one_lt_rpow_iff_of_pos (lt_trans zero_lt_one h1) |>.mpr
    (Or.inl ⟨h1, div_pos one_pos (by linarith)⟩)

/-- ★★★**深さ 1 だけ `axDecay p 1` に落としても、有限積は非有界**。

⇒ ★★`AxLemma K C` はどんな `C` についても出ない
（`axLemma_of_wildDescent_Icc` の仮説 `∀ n, ∏_{Icc 1 n} c k ≤ C` が満たせない）。
★**すべての深さ `k` で `axDecay p k` が要る**ことの証明である。

★`PureStepSetup` の道は `k = 0` 専用（同ファイル `:430`）なので、
★★あの構造を 1 つ作っても届くのは深さ 1 だけで、そこで止まる。 -/
theorem prod_depth_one_only_unbounded (C : ℝ) :
    ∃ n : ℕ, C < ∏ k ∈ Finset.Icc 1 n,
      (if k = 1 then axDecay p 1 else (p : ℝ)) := by
  have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have ha : (1 : ℝ) < axDecay p 1 := one_lt_axDecay_one
  refine prod_unbounded_of_one_lt (D := min (axDecay p 1) (p : ℝ)) (lt_min ha h1) ?_ C
  intro k
  by_cases hk : k = 1
  · simp [hk]
  · simp [hk]

end NotEnough

/-! ## §3 使っている公理の一覧 -/

#print axioms axWildDescent_of_wildStep_approx
#print axioms one_lt_axDecay_one
#print axioms prod_depth_one_only_unbounded

end WildStepFieldSupply

end ABC3.Found.PGC
