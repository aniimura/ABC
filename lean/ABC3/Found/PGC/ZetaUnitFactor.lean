import ABC3.Found.PGC.ZetaStepRatio

/-!
# [pGC] ★★★★柱が**仮定ゼロ**になった ＋ ★★**棚卸しの更新**（次の波はここを最初に読む）

## ★どれを選んだか、なぜか（費用を先に測った）

候補 (P)「`‖·‖ ↔ v_L` の記法翻訳」を測った。★**測った結果、翻訳は要らなかった。**

```
grep -rn "Valuation|valuation" （本日の 7 ファイル）
  ⇒ Lean の**文**に `Valuation` は 1 つも無い。全部 `‖·‖` で書かれている。
    当たったのは私自身が置いた**数値の置き石** 2 つだけ:
      RhoFactorization.lean:174 w_valuation_p      : (2:ℕ)=2 ∧ (3:ℕ)=3 ∧ (5:ℕ)=5
      RhoFactorization.lean:179 sigma_mu_valuation_sq : 4=2^2 ∧ 9=3^2 ∧ 25=5^2
```

⇒ ★★★**訂正: 前波で私が「残るのは `v_L` への翻訳（配管）」と書いたのは誤り。**
★実際に残っていたのは、上の 2 つが**まだ測定のまま**だったことである。
（★前波は「1 行」の見積もりを外し、今回は「翻訳が要る」を外した。★2 度目である。）

⇒ 本波はその 2 つを**定理に上げた**。所要 2 往復。
★もう 1 つの候補 (Q)「`‖p‖ < 1` の instance」は #69 の境界に近いので今回も選ばず、
(A) ①(a) は測定済みで重い、(G)/(R) は本波の流れから遠い、と判断した。
★(S)（棚卸しの更新）は**同じ波で一緒にやった**（下の §0）。

## §0 ★★★棚卸し（2026-09-09 現在。`PairBudgetVerified.lean` の表は **12 波ぶん古い**）

### 穴の現状

| 穴 | 状態 | 根拠 |
|---|---|---|
| ①不分岐 | (b) **閉**／(a) 残（★重いと測定済み） | `UnramifiedStepFixed` / `UnramifiedLayerFree` |
| ②′一様定数 | ★**測れた** —— 目盛りで測ると定数は深さに依らない | `HformPopulation` / `LossUnitLaw` |
| ③幾何減衰 | ★**支持された** —— 鋭い定数も公比 `1/p` の等比列 | `HformPopulation.sharp_geometric` |
| ④測定点 | **消えた**（積では通る） | `BudgetFiniteExcess` |
| ⑤出口の限界 | ★**確定** —— `axDecay` は 1 目盛り足りない | `LossUnitLaw` / `FirstJumpNotAchieved` |

### `hform` の鎖（12 波。★ここが本日の主線）

| # | ファイル | 結論 |
|---|---|---|
| 1 | `FirstJumpNotAchieved` | `hform` は測定点で偽。表を初再現、Herbrand 同定 `i₁ = j = 2` |
| 2 | `HformPopulation` | 99.7% で成立。★核 `axDecay p k = p^{p/φ(p^{k+1})}`（層の付値で**ちょうど `p` 目盛り**） |
| 3 | `LossUnitLaw` | ★配られた `loss ≤ 2·i₁` は **`p=2` で偽**。深さ 1 は安全（`H_E = H_K`） |
| 4 | `FirstJumpWitness` | 証人 1 個への還元。★最上位成分だけの仮説は**自分で潰した** |
| 5 | `TailNoCancel` | ★★**打ち消しの原因は `f₀`** と確定（尾だけ 16,500 件で分散ゼロ）。`hform` を**初めて導いた**。必要条件 `d ≡ 1 (mod p)` を**証明** |
| 6 | `SlotStructure` | ★深さ ＝ 生き残る成分の番号 `j`（空き枠の補題を**証明**） |
| 7 | `LossPAddOneFalse` | ★★`loss ≤ p+1` は **`p=5` で偽**（明示的反例、7 秒）。`p+1` は**2 つの別々の偶然** |
| 8 | `TwoIsDegenerate` | `p=2` が特別な理由 2 つを**証明**（`p(p−1)=p ⟺ p=2`、尾が 1 本） |
| 9 | `RhoFactorization` | ★2 本の柱が**同じ 1 本の恒等式** `ρ = (1+π)·w` に落ちた |
| 10 | `ZetaSubOnePrime` | ★最後の仮定「`ζ−1` が素元」が**定理**に（mathlib Eisenstein ＋ 木の `norm_pow_eq_of_monic_root`） |
| 11 | `ZetaStepRatio` | ★`v_L(μ) = p` が**定理** |
| 12 | ★本ファイル | ★★`v_L(w) = p` と `v_L(σμ−μ) = p²` が**定理** ⇒ **柱は仮定ゼロ** |

### 掘らなくてよいと確定した道（4 本）

`AxWildDescentDecay`（`AxLemma` と同値＝循環）／予算関数 `F` への載せ替え（同じ循環）／
`FirstJumpRoute` の再パラメータ化（`c k ≤ axDecay p k` と同値）／B2（発火しない）。

### ★本日の「小さい場合の偶然の一致」（#350）—— **4 件**

`zeta16_sharp` の目盛り／`4 = 2·i₁`／`loss ≤ p+1`（★2 つの偶然が重なった）／`t = e_{E₁}`。
⇒ ★**法則を立てたら必ず 3 つ目の場合で測る。**

## §1〜§4 の内容（本波の新規）

- `norm_root_of_unity_one` —— 単位根のノルムは 1。
- `norm_pow_sub_one_le` —— ★`‖ξ^i − 1‖ ≤ ‖ξ − 1‖`（帰納法。`ξ^{i+1}−1 = ξ(ξ^i−1)+(ξ−1)` と強三角不等式だけ）。
- `norm_geom_sum_sub_natCast_le` —— `‖Σ_{i<b} ξ^i − b‖ ≤ ‖ξ − 1‖`。
- `norm_geom_sum_eq_one` —— ★**幾何和は単数**（`‖(b:F)‖ = 1` すなわち `p ∤ b` のとき）。
- `norm_pow_sub_one_eq` —— ★★**`‖ξ^b − 1‖ = ‖ξ − 1‖`** ＝ `v_L(w) = v_L(μ)`。
  ★これが `RhoFactorization` の `w_valuation_p`（測定）を**定理にした**もの。
- `norm_zeta_step_two` —— ★★`(‖ζ−1‖^p)^p = ‖ξ−1‖` ＝ `v_L(σμ−μ) = p²`。
  ★これが `sigma_mu_valuation_sq`（測定）を**定理にした**もの。

## ★残り（★正直に）

1. ★`‖p‖ < 1` と `‖(b:F)‖ = 1`（`p ∤ b`）は**仮定のまま**である。
   どちらも「`F` の剰余標数が `p`」から出るが、その instance は与えていない（#69 の境界）。
2. ★①(a)、②′の証人、③、⑤ は**変化なし**。
3. ★`GainedTowerDescent.lean:94` の総当たりの表（掃討の残り 1 箇所）は**未着手**。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
   ★`RhoFactorization` の 2 つの置き石は**そのまま残す** —— 他が読んでいるので直さず、
   ★**本ファイルが上位互換であることをここに書く**（CLAUDE.md「訂正は自分のファイルに」）。
3. ★本波は **(S) 棚卸しの更新を同じ波で兼ねた**（本日この鎖は圧縮を 3 回挟んでおり、
   ★圧縮後の自分が最初に読む場所を最新にする価値が高いと判断した）。
4. ★宣言名 6 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
5. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace ZetaUnitFactor

/-! ## §1 単位根のノルムは 1 -/

section RootOfUnity

/-- ★単位根のノルムは 1（`‖ζ‖^N = ‖ζ^N‖ = 1` と冪の消去）。 -/
theorem norm_root_of_unity_one {F : Type*} [NormedField F] {ζ : F} {N : ℕ}
    (hN : 0 < N) (h : ζ ^ N = 1) : ‖ζ‖ = 1 := by
  refine ZetaStepRatio.pow_left_cancel_nonneg (norm_nonneg _) zero_le_one hN ?_
  rw [one_pow, ← norm_pow, h, norm_one]

end RootOfUnity

/-! ## §2 `‖ξ^i − 1‖ ≤ ‖ξ − 1‖` -/

section PowSubOne

/-- ★`‖ξ^i − 1‖ ≤ ‖ξ − 1‖`。
★帰納法 —— `ξ^{i+1} − 1 = ξ(ξ^i − 1) + (ξ − 1)` と強三角不等式だけ。
★幾何和の因数分解を経由するより短い。 -/
theorem norm_pow_sub_one_le {F : Type*} [NormedField F] [IsUltrametricDist F]
    {ξ : F} (hξ : ‖ξ‖ ≤ 1) : ∀ i : ℕ, ‖ξ ^ i - 1‖ ≤ ‖ξ - 1‖
  | 0 => by simp
  | (i + 1) => by
      have hstep : ξ ^ (i + 1) - 1 = ξ * (ξ ^ i - 1) + (ξ - 1) := by ring
      have hih := norm_pow_sub_one_le hξ i
      calc ‖ξ ^ (i + 1) - 1‖
          = ‖ξ * (ξ ^ i - 1) + (ξ - 1)‖ := by rw [hstep]
        _ ≤ max ‖ξ * (ξ ^ i - 1)‖ ‖ξ - 1‖ := IsUltrametricDist.norm_add_le_max _ _
        _ ≤ ‖ξ - 1‖ := by
            refine max_le ?_ le_rfl
            rw [norm_mul]
            calc ‖ξ‖ * ‖ξ ^ i - 1‖ ≤ 1 * ‖ξ ^ i - 1‖ := by
                  exact mul_le_mul_of_nonneg_right hξ (norm_nonneg _)
              _ = ‖ξ ^ i - 1‖ := one_mul _
              _ ≤ ‖ξ - 1‖ := hih

end PowSubOne

/-! ## §3 幾何和は単数 -/

section GeomSum

/-- `‖Σ_{i<b} ξ^i − b‖ ≤ ‖ξ − 1‖`（`Σ ξ^i − b = Σ (ξ^i − 1)` と超距離の和）。 -/
theorem norm_geom_sum_sub_natCast_le {F : Type*} [NormedField F] [IsUltrametricDist F]
    {ξ : F} (hξ : ‖ξ‖ ≤ 1) (b : ℕ) :
    ‖(∑ i ∈ Finset.range b, ξ ^ i) - (b : F)‖ ≤ ‖ξ - 1‖ := by
  have hrw : (∑ i ∈ Finset.range b, ξ ^ i) - (b : F)
      = ∑ i ∈ Finset.range b, (ξ ^ i - 1) := by
    rw [Finset.sum_sub_distrib]
    simp
  rw [hrw]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (norm_nonneg _) ?_
  intro i _
  exact norm_pow_sub_one_le hξ i

/-- ★★**幾何和は単数** —— `‖(b:F)‖ = 1`（＝ `p ∤ b`）かつ `‖ξ − 1‖ < 1` なら
`‖Σ_{i<b} ξ^i‖ = 1`。★強三角不等式の**等号版**（`FirstJumpWitness.norm_add_eq_of_norm_lt`）。 -/
theorem norm_geom_sum_eq_one {F : Type*} [NormedField F] [IsUltrametricDist F]
    {ξ : F} {b : ℕ} (hξ : ‖ξ‖ ≤ 1) (hb : ‖(b : F)‖ = 1) (hlt : ‖ξ - 1‖ < 1) :
    ‖∑ i ∈ Finset.range b, ξ ^ i‖ = 1 := by
  have hkey : ‖((∑ i ∈ Finset.range b, ξ ^ i) - (b : F)) + (b : F)‖ = ‖(b : F)‖ :=
    FirstJumpWitness.norm_add_eq_of_norm_lt
      (lt_of_le_of_lt (norm_geom_sum_sub_natCast_le hξ b) (by rw [hb]; exact hlt))
  rw [sub_add_cancel, hb] at hkey
  exact hkey

/-- ★★★**`‖ξ^b − 1‖ = ‖ξ − 1‖`** ＝ `v_L(w) = v_L(μ)`。

★`RhoFactorization.w_valuation_p`（3 素数での**測定**）を**定理にした**もの。
⇒ 本日の柱の 1 本目（`ρ = (1+π)·w` は 2 成分）が仮定ゼロになる。 -/
theorem norm_pow_sub_one_eq {F : Type*} [NormedField F] [IsUltrametricDist F]
    {ξ : F} {b : ℕ} (hξ : ‖ξ‖ ≤ 1) (hb : ‖(b : F)‖ = 1) (hlt : ‖ξ - 1‖ < 1) :
    ‖ξ ^ b - 1‖ = ‖ξ - 1‖ := by
  have hfac : ξ ^ b - 1 = (ξ - 1) * ∑ i ∈ Finset.range b, ξ ^ i :=
    RhoFactorization.geom_unit_factor ξ b
  rw [hfac, norm_mul, norm_geom_sum_eq_one hξ hb hlt, mul_one]

end GeomSum

/-! ## §4 2 段の比 -/

section TwoStep

/-- ★★★`(‖ζ − 1‖^p)^p = ‖ξ − 1‖` ＝ `v_L(σμ − μ) = p²`。

★`RhoFactorization.sigma_mu_valuation_sq`（3 素数での**測定**）を**定理にした**もの。
⇒ 本日の柱の 2 本目（`t = p(p−1)`）が仮定ゼロになる。 -/
theorem norm_zeta_step_two {F : Type*} [NormedField F] [IsUltrametricDist F]
    [CharZero F] {p m : ℕ} [Fact p.Prime] {ζ η ξ : F}
    (hζ : IsPrimitiveRoot ζ (p ^ (m + 3))) (hη : IsPrimitiveRoot η (p ^ (m + 2)))
    (hξ : IsPrimitiveRoot ξ (p ^ (m + 1)))
    (hlt : ‖(p : F)‖ < 1) (hne : (p : F) ≠ 0) :
    (‖ζ - 1‖ ^ p) ^ p = ‖ξ - 1‖ := by
  rw [ZetaStepRatio.norm_zeta_step (m := m + 1) hζ hη hlt hne]
  exact ZetaStepRatio.norm_zeta_step (m := m) hη hξ hlt hne

end TwoStep

/-! ## §5 使っている公理の一覧 -/

#print axioms norm_root_of_unity_one
#print axioms norm_pow_sub_one_le
#print axioms norm_geom_sum_sub_natCast_le
#print axioms norm_geom_sum_eq_one
#print axioms norm_pow_sub_one_eq
#print axioms norm_zeta_step_two

end ZetaUnitFactor

end ABC3.Found.PGC
