import ABC3.Found.PGC.WildJumpChain
import ABC3.Found.PGC.RamificationJumpBound

/-!
# [pGC] 台帳の検算 —— 鎖による降下が `axDecay p k` に届くのは「最後の跳びが最大」のときだけ

`Found/PGC/WildJumpChain.lean` は `AxWildDescent K (axDecay p)` に残るものを 3 点
(1) `(p−1)u₁ ≤ p·e_K`、(2) `hbase`(`m = 1` を収縮率の言葉に翻訳)、
(3) 最初の跳びの層で降ろしても wild 深さが真に下がること、と名指ししていた。

★★★本ファイルの結論は 4 つである。★★埋まっていないものは埋まっていないと書く。

* ★★★**台帳の勘定は閉じない**。しかも「閉じない例が在る」より強く
  ★**閉じるのは `s_m = 1`(最後の跳びが最大 `(p−1)i = e_L`)のとき★だけ★**である
  (`chain_ledger_forces_max`、★実数だけ、★`sorry` 無し)。
  ★系として `p ∣ i` が強制される(`dvd_of_chain_ledger`)ので、
  ★`p ∤ i` の実在データでは★端から閉じない★(`not_chain_ledger_of_not_dvd`)。
* ★★点 2 は★木に既に在った★。`RamificationJumpBound.norm_natCast_le_pow_of_splits` の
  結論 `‖(n:L)‖ ≤ ‖π‖^{(n−1)i}` が `n = p` でそのまま `‖p‖ ≤ (‖π‖^i)^{p−1}` である
  (橋は `pow_mul` の並べ替え 1 行 = `norm_natCast_le_pow_of_mul_le`)。
  ★★ただし `WildJumpChain` の `hbase` は `σ : K.absGal` と `∀ z : K.closure` で書かれており、
  ★その形は下の §測定 3 のとおり**満たせない**。
* ★点 3 の前提が外れている。★**「ちょうど 1 下がる」を要求している箇所は無い**
  (`AxWildDescent` も `AxTowerDecay.exists_mem_of_descent_budget` も `<` しか要求しない)。
  ★しかも降下自体は `WildDepthFieldDescent.axWildDescent_normInv` が**既に持っている**。
  ⇒ ★配置替えは要らない。足りないのは★定数だけ★である(`axWildDescent_mono`)。
* ★点 1 は★閉じなかった★。落ちたのは `ℕ` の糊(`sub_one_mul_first_jump_le`)と
  ★測って分かった弱い版(`sub_one_mul_first_jump_le_of_metric`)だけで、
  sharp な `(p−1)u₁ ≤ p·e_K` は Herbrand が要る(§測定 5)。

## ★★★測定 1 —— 台帳の正しい字面(★`WildJumpChain` の字面は損失を過小評価している)

`σ_j = τ^{p^j}`(`j = 0..m`)の収縮率を `θ_j`、跳びを `i_j`、`e_L = v_L(p)` と置き

  `s_j := (p−1)·i_j / e_L`   (`θ_j = ‖π_L‖^{i_j}`、`‖p‖ = ‖π_L‖^{e_L}`)

とする。`CyclicLayerDescent` が実際に持っている 2 本を指数に直すと

| 補題 | 値 | 指数(`p` 進、`p^{−·}`) |
|---|---|---|
| `norm_sub_orbit_average_le_of_contract` | 損失 `max (p·θ_m^{p−2}) 1` | ★`1 − (p−2)s_m/(p−1)` |
| `norm_smul_pow_prime_sub_self_le_of_contract` | 縮み `max (θ_j^{p−1}) ‖p‖` | `−min(s_j, 1) = −s_j` |

⇒ 総損失の指数は ★`1 − (p−2)·s_m/(p−1) − Σ_{j<m} s_j`、
`axDecay p (m+1)` の指数は `p^{−m}/(p−1)`。これが本ファイルの `hledger` の字面である。

★★`WildJumpChain` の `ledger_le` / `not_ledger_of_ge` は損失を `s_m/(p−1)` と書いていた。
★2 つが一致するのは `s_m = 1` のとき★だけ★で、`s_m < 1` では
`1 − (p−2)s_m/(p−1) > s_m/(p−1)`、すなわち★真の損失はもっと大きい★。
(実在データ `ℚ₃(ζ₂₇)` で `WildJumpChain` は `2/9 > 1/6` と測ったが、正しくは `1/3 > 1/6`。)

## ★★★★測定 2 —— 剛性: 閉じるなら `s_m = 1`、したがって `p ∣ i_m`

鎖から出るのは `s_{j+1} ≥ p·s_j`(★下から。`WildJumpChain` §測定 3)なので
`s_j ≤ p^{j−m}·s_m`、よって `Σ_{j<m} s_j ≤ s_m·(1−p^{−m})/(p−1)`
(`sum_le_of_chain_ge`)。これを `hledger` に入れて `(p−1)` を払うと

  `(p−1) − p^{−m} ≤ s_m·((p−2) + 1 − p^{−m}) = s_m·(p−1−p^{−m})`

で `p−1−p^{−m} > 0`(`m ≥ 1`)だから ★`s_m ≥ 1`。`s_m ≤ 1`(= `RamificationJumpBound`)と
合わせて ★★`s_m = 1`、すなわち `(p−1)·i_m = e_L = p·e_F`(`F` は `σ_m` の固定体)。

`RamificationJumpBound.sub_one_mul_lt_of_not_dvd` は「`p ∤ i` なら `(p−1)i < p·e`」なので

  ★★★勘定が閉じる ⟹ `p ∣ i_m`  (`dvd_of_chain_ledger`)

★実在データ `p = 3`、`L = ℚ₃(ζ₂₇)`、`K = ℚ₃(ζ₃)`、`F = ℚ₃(ζ₉)`:
`i(τ³) = 8`、`e_L = 18 = 3·6 = p·e_F`、`s_m = 16/18 = 8/9 < 1`、★`3 ∤ 8`。
⇒ ★★このデータでは端から閉じない(`not_chain_ledger_of_not_dvd`、数値版は
`not_chain_ledger_cyclotomic` で `1/3 ≤ 1/6` が偽であることを示した)。

★★つまり「鎖の勘定」は Kummer 塔(`K = ℚ_p(ζ_{p²})`, `x = p^{1/p²}`、跳びが等比)という
★特殊な 1 系列でしか閉じない。★★★**3 点が全部埋まっても `axDecay p k` は出ない。**

## ★★★★測定 3 —— 収縮率の仮定そのものが空虚である(★★形式化していない。手計算)

`CyclicLayerDescent` / `WildJumpChain` の `_of_contract` 族は

  `hθ : ∀ z : K.closure, ‖σ • z − z‖ ≤ θ * ‖z‖`

を仮定する。★★これは `θ < 1` では★満たせない★。逐語の証人(手計算):

  `K = ℚ_p`、`σ ∈ G_K` を円分指標が `χ(σ) = 1 + p^n`(`n ≥ 1`)となるものに取る。
  `z = ζ_{p^k} − 1`(`k > n`)に対し `σz − z = ζ_{p^k}(ζ_{p^{k−n}} − 1)` なので
  `‖σz − z‖ / ‖z‖ = p^{−(p^n − 1)/(p^{k−1}(p−1))} → 1`  (`k → ∞`)。

⇒ ★`sup_z ‖σz − z‖/‖z‖ = 1`。★一般に `K̄/K` は deeply ramified(Coates–Greenberg)なので
`σ ≠ 1` なら常に `1` である。★★したがって `hθ` を満たすのは `θ = 1`(または `σ = 1`)だけで、
そこでは道が★退化する★:

  損失 `max (p·1^{p−2}) 1 = p`(`orbit_average_loss_at_one`)、
  縮み `max (1^{p−1}) ‖p‖ = 1`(`contract_gain_at_one`)

—— ★`WildDepthFieldDescent.axWildDescent_prime` の定数 `p` そのものに戻る。

★★★同じ理由で `WildJumpChain` の `hbase` も字面のままでは使えない。加えて
★**`K.absGal` に位数 `p` の元は無い**(Artin–Schreier: 絶対 Galois 群の非自明な有限部分群は
位数 2 で固定体が実閉体。`p` 進体は形式的実でない)ので、
`hbase` の根拠「`σ^{p^m}` は位数 `p`」は `K.absGal` では★決して成り立たない★。
★`σ = 1` を入れると `η = 0` が取れて `‖p‖ ≤ 0` になるから `hbase` は★偽★でもある。

★★★直し方(測ったので書いておく): `_of_contract` の証明は `D = σ−1` の反復しか使わず、
`z` は `K(x)` の Galois 閉包という★有限層★から出ない。
⇒ 仮定を「`∀ z ∈ M`(有限層)」に弱めれば `θ = ‖π_M‖^{i_M(σ)} < 1` が正当になる。
★本ファイルはその書き換えを行っていない(持ち場が 4 本を「読むだけ」としたため)。

## ★測定 4 —— 点 3: 「ちょうど 1 下がる」を要求している箇所は無い

`AxWildDescent K c` の結論は `wildDepth K x' < wildDepth K x` で★狭義の減少だけ★、
`AxTowerDecay.exists_mem_of_descent_budget` の再帰も `hlt : deg x' < n` で強い帰納法に入る。
`exists_mem_of_descent_prod` の予算 `F d = ∏_{k ∈ Icc 1 d} c k` は `1 ≤ c k` があるので
★段を飛ばしても損しない★。⇒ ★`(ℤ/p)^k` で一気に深さ 0 まで落ちても壊れない。
★さらに降下自体は `axWildDescent_normInv` が既に(深さちょうど 1 減で)持っている。
⇒ ★★点 3 は「新規に要る点」ではない。定数だけが問題である(`axWildDescent_mono`)。

## ★測定 5 —— 点 1 は閉じなかった(★理由を正確に書く)

要るのは `u₁`(`L/K` の最初の跳び)について ★`(p−1)·u₁ ≤ p·e_K`(★底 `K`)。

* `ℕ` の糊は落ちる: 指数 `p` の商 `L^N/K` の跳び `i` について `u₁ ≤ i` と
  `(p−1)i ≤ p·e_K` を繋ぐだけ(`sub_one_mul_first_jump_le`)。
* ★足りないのは `u₁ ≤ i(σ̄)`(`σ̄` は `G/N` の生成元)である。これは
  Herbrand(`φ_{L/K}` が `u ≤ u₁` で恒等)、同値に Serre, Corps Locaux IV §1 Prop. 3 の
  `i_{L^N}(σ̄) = |N|^{−1}·Σ_{ν∈N} i_L(σν)` である。★本ファイルは形式化していない。
* ★★ノルムだけで代用しようとすると★弱くなる★ことを測った: `‖σ̄w − w‖ = ‖σw − w‖` から
  出るのは `i(σ̄)/e_{L^N} ≥ u₁/e_L`(正規化された比)であって `i(σ̄) ≥ u₁` ではない。
  そこから出る結論は ★`(p−1)u₁ ≤ e_L`(`sub_one_mul_first_jump_le_of_metric`)で、
  ★欲しい `(p−1)u₁ ≤ p·e_K` より `[L:K]/p` 倍だけ弱い。
  ⇒ ★「距離だけでは Herbrand を代用できない」ことが測れた。

★なお ★測定 2 のとおり台帳が閉じないので、点 1 を閉じても `axDecay p k` は出ない。

## 在庫の測定(★自分で測った。★MCP は 1 度も呼んでいない。★引いたコマンドを残す)

```
grep -c "\b<名前>\b" .cache/decl-index.txt   （本ファイルの 12 宣言すべて 0 = 衝突なし）
sed -n '316,345p' lean/ABC3/Found/PGC/RamificationJumpBound.lean
  → ★★norm_natCast_le_pow_of_splits の結論は ‖(n : L)‖ ≤ ‖π‖ ^ ((n - 1) * i)。
     ★これは n = p でそのまま hbase(‖p‖ ≤ (‖π‖^i)^{p−1})である。
     ★「木に無い」と WildJumpChain が書いた点 2 は★在った★——
     ★探し方が「収縮率 θ の名前」で、実際は「π と i の形」で持っていたためである。
grep -in "deeplyramified|deeply_ramified" .cache/mathlib-index.txt   → 0 件
grep -in "ArtinSchreier|artin_schreier"   .cache/mathlib-index.txt   → 0 件
  ⇒ ★測定 3 を機械化する語彙は mathlib に無い(測ったうえで「無い」と書く)。
```

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★配られた 3 点のうち閉じたのは点 2(既に木に在った)だけである。点 1 は閉じず、
   点 3 は★そもそも要らない★ことが測れた。★どれも報告のとおりには落ちていない。
2. ★★配られた「3 点で足りる」は★偽★である(測定 2)。3 点が全部埋まっても
   鎖の勘定は `s_m = 1` 以外で閉じない。★これが本ファイルの主結果である。
3. ★`WildJumpChain` の `ledger_le` / `not_ledger_of_ge` の損失項 `s_m/(p−1)` は
   `s_m < 1` で過小評価である(測定 1)。★あちらは真の主張だが、
   ★台帳の字面としては本ファイルの `hledger` の方が正しい。★あちらは書き換えていない。
4. ★測定 3(収縮率の仮定の空虚性)は形式化していない。手計算であることを明記した。
5. `.src` は `CyclicLayerDescent.lean` / `WildJumpChain.lean` と同じ項目
   (pGC 物理 p.6 Corollary 3.1)を指す。原典が独立に立てた項目ではない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

/-! ## §1 抽象核(実数だけ) -/

section Chain

/-- 抽象核(実数だけ): `p·s_j ≤ s_{j+1}` を `d` 段まわす。

★`WildJumpChain.step_le`(`s_{j+1} ≤ p·s_j`)の★向き違い★である。
鎖(`norm_smul_pow_prime_sub_self_le_of_contract`)から実際に出るのはこちらの向き。 -/
theorem step_ge {p : ℕ} {s : ℕ → ℝ} {m : ℕ}
    (hrec : ∀ j, j < m → (p : ℝ) * s j ≤ s (j + 1)) (hp0 : (0:ℝ) ≤ (p:ℝ)) :
    ∀ (d j : ℕ), j + d ≤ m → (p : ℝ) ^ d * s j ≤ s (j + d) := by
  intro d
  induction d with
  | zero => intro j _; simp
  | succ e ih =>
      intro j hj
      have h2 : (p:ℝ) ^ e * s j ≤ s (j + e) := ih j (by omega)
      have h1 : (p:ℝ) * s (j + e) ≤ s (j + e + 1) := hrec (j + e) (by omega)
      have hje : (j + (e + 1)) = (j + e + 1) := by omega
      rw [hje]
      calc (p:ℝ) ^ (e + 1) * s j = (p:ℝ) * ((p:ℝ) ^ e * s j) := by ring
        _ ≤ (p:ℝ) * s (j + e) := mul_le_mul_of_nonneg_left h2 hp0
        _ ≤ s (j + e + 1) := h1

def step_ge.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- 抽象核(実数だけ): 鎖 `p·s_j ≤ s_{j+1}` は総和の★上界★を与える。

  `Σ_{j<m} s_j ≤ ((1 − p^{−m})/(p−1))·s_m`

★これが「鎖から出る情報は稼ぎの★上限★でしかない」ことの本体である。
台帳は稼ぎの★下限★を要求するので、ここで向きが噛み合わない。 -/
theorem sum_le_of_chain_ge {p : ℕ} (hp : 1 < (p : ℝ)) {s : ℕ → ℝ} (m : ℕ)
    (hrec : ∀ j, j < m → (p : ℝ) * s j ≤ s (j + 1)) :
    (∑ j ∈ Finset.range m, s j) ≤ ((1 - ((p:ℝ)⁻¹) ^ m) / ((p : ℝ) - 1)) * s m := by
  have hp0 : (0:ℝ) < (p : ℝ) := lt_trans zero_lt_one hp
  have hkey : ∀ j ∈ Finset.range m, s j ≤ ((p:ℝ)⁻¹) ^ (m - j) * s m := by
    intro j hj
    rw [Finset.mem_range] at hj
    have h := step_ge hrec (le_of_lt hp0) (m - j) j (by omega)
    have hmj : j + (m - j) = m := by omega
    rw [hmj] at h
    rw [inv_pow, ← div_eq_inv_mul, le_div_iff₀ (by positivity)]
    calc s j * (p:ℝ) ^ (m - j) = (p:ℝ) ^ (m - j) * s j := by ring
      _ ≤ s m := h
  have hsum : (∑ j ∈ Finset.range m, s j)
      ≤ (∑ j ∈ Finset.range m, ((p:ℝ)⁻¹) ^ (m - j)) * s m := by
    rw [Finset.sum_mul]
    exact Finset.sum_le_sum hkey
  have hgeom : (∑ j ∈ Finset.range m, ((p:ℝ)⁻¹) ^ (m - j))
      = (1 - ((p:ℝ)⁻¹) ^ m) / ((p:ℝ) - 1) := by
    rw [← sum_inv_pow_succ_eq hp m,
      ← Finset.sum_range_reflect (fun i => ((p:ℝ)⁻¹) ^ (i + 1)) m]
    refine Finset.sum_congr rfl (fun j hj => ?_)
    rw [Finset.mem_range] at hj
    congr 1
    omega
  rwa [hgeom] at hsum

def sum_le_of_chain_ge.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★★**剛性**(本ファイルの主結果。★実数だけ、★分岐も付値も Galois も出てこない)。

鎖の仮定(`p·s_j ≤ s_{j+1}`、`s_m ≤ 1`)のもとで台帳

  `1 − (p−2)·s_m/(p−1) − Σ_{j<m} s_j ≤ p^{−m}/(p−1)`

が成り立つのは ★`s_m = 1` のとき★だけ★である。

★★つまり「鎖による勘定」は★最後の跳びが最大 `(p−1)i_m = e_L` の場合しか閉じない★。
`m ≥ 1`(= 深さ `k ≥ 2`)が要る: `m = 0` は `axDecay p 1` の定義そのもので中身が無い。 -/
theorem chain_ledger_forces_max {p : ℕ} (hp : 1 < (p : ℝ)) {s : ℕ → ℝ} {m : ℕ} (hm : 1 ≤ m)
    (hsm : s m ≤ 1)
    (hrec : ∀ j, j < m → (p : ℝ) * s j ≤ s (j + 1))
    (hledger : 1 - ((p : ℝ) - 2) * s m / ((p : ℝ) - 1) - (∑ j ∈ Finset.range m, s j)
        ≤ ((p : ℝ)⁻¹) ^ m / ((p : ℝ) - 1)) :
    s m = 1 := by
  have hp0 : (0:ℝ) < (p:ℝ) := lt_trans zero_lt_one hp
  have hpn : 1 < p := by exact_mod_cast hp
  have hp2 : (2:ℝ) ≤ (p:ℝ) := by exact_mod_cast hpn
  have hP1 : (0:ℝ) < (p:ℝ) - 1 := by linarith
  have hP1ne : ((p:ℝ) - 1) ≠ 0 := ne_of_gt hP1
  have hq0 : (0:ℝ) < ((p:ℝ)⁻¹) := by positivity
  have hinv : ((p:ℝ)⁻¹) * (p:ℝ) = 1 := inv_mul_cancel₀ (ne_of_gt hp0)
  have hq1 : ((p:ℝ)⁻¹) ≤ 1/2 := by nlinarith [hq0, hp2, hinv]
  have hqm : ((p:ℝ)⁻¹) ^ m ≤ ((p:ℝ)⁻¹) :=
    pow_le_of_le_one (le_of_lt hq0) (by linarith) (by omega)
  have hqm2 : ((p:ℝ)⁻¹) ^ m ≤ 1/2 := le_trans hqm hq1
  have hsum := sum_le_of_chain_ge hp m hrec
  have h1 : ((p:ℝ) - 1) - ((p:ℝ) - 2) * s m
      - ((p:ℝ) - 1) * (∑ j ∈ Finset.range m, s j) ≤ ((p:ℝ)⁻¹) ^ m := by
    have hmul := mul_le_mul_of_nonneg_left hledger (le_of_lt hP1)
    have e1 : ((p:ℝ) - 1) * (1 - ((p:ℝ) - 2) * s m / ((p:ℝ) - 1)
          - (∑ j ∈ Finset.range m, s j))
        = ((p:ℝ) - 1) - ((p:ℝ) - 2) * s m
          - ((p:ℝ) - 1) * (∑ j ∈ Finset.range m, s j) := by
      field_simp
    have e2 : ((p:ℝ) - 1) * (((p:ℝ)⁻¹) ^ m / ((p:ℝ) - 1)) = ((p:ℝ)⁻¹) ^ m := by
      field_simp
    rw [e1, e2] at hmul
    exact hmul
  have h2 : ((p:ℝ) - 1) * (∑ j ∈ Finset.range m, s j)
      ≤ (1 - ((p:ℝ)⁻¹) ^ m) * s m := by
    have hmul := mul_le_mul_of_nonneg_left hsum (le_of_lt hP1)
    have e3 : ((p:ℝ) - 1) * (((1 - ((p:ℝ)⁻¹) ^ m) / ((p:ℝ) - 1)) * s m)
        = (1 - ((p:ℝ)⁻¹) ^ m) * s m := by
      field_simp
    rw [e3] at hmul
    exact hmul
  have hP1q : (0:ℝ) < (p:ℝ) - 1 - ((p:ℝ)⁻¹) ^ m := by linarith
  have hge : 1 ≤ s m := by nlinarith [h1, h2, hP1q]
  exact le_antisymm hsm hge

def chain_ledger_forces_max.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Chain

/-! ## §2 付値語への翻訳 —— 勘定が閉じるなら `p ∣ i` -/

/-- ★★★剛性の付値語版: `s_m = (p−1)i/(p·e_F)` を入れると ★`p ∣ i` が★強制される★。

`RamificationJumpBound.sub_one_mul_lt_of_not_dvd`(`p ∤ i` なら `(p−1)i < p·e`)の対偶。
★`e_L = p·e_F`(`F` は位数 `p` の元 `σ_m` の固定体)である。 -/
theorem dvd_of_chain_ledger {p i eK : ℕ} (hp : 1 < p) (heK : 0 < eK) {s : ℕ → ℝ} {m : ℕ}
    (hm : 1 ≤ m)
    (hs : s m = (((p:ℝ) - 1) * (i : ℝ)) / ((p : ℝ) * (eK : ℝ)))
    (hsm : s m ≤ 1)
    (hrec : ∀ j, j < m → (p : ℝ) * s j ≤ s (j + 1))
    (hledger : 1 - ((p : ℝ) - 2) * s m / ((p : ℝ) - 1) - (∑ j ∈ Finset.range m, s j)
        ≤ ((p : ℝ)⁻¹) ^ m / ((p : ℝ) - 1)) :
    p ∣ i := by
  have hpR : (1:ℝ) < (p:ℝ) := by exact_mod_cast hp
  have hp0 : (0:ℝ) < (p:ℝ) := lt_trans zero_lt_one hpR
  have heK0 : (0:ℝ) < (eK:ℝ) := by exact_mod_cast heK
  have hone := chain_ledger_forces_max hpR hm hsm hrec hledger
  rw [hs, div_eq_one_iff_eq (by positivity)] at hone
  have hcast : (((p - 1) * i : ℕ) : ℝ) = ((p * eK : ℕ) : ℝ) := by
    push_cast [Nat.cast_sub (le_of_lt hp)]
    exact hone
  have hnat : (p - 1) * i = p * eK := by exact_mod_cast hcast
  by_contra hnd
  have hlt := sub_one_mul_lt_of_not_dvd (by omega : 0 < p) hnd (le_of_eq hnat)
  omega

def dvd_of_chain_ledger.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★対偶: `p ∤ i_m` なら鎖の勘定は★端から閉じない★。

★実在データ: `p = 3`、`L = ℚ₃(ζ₂₇)`、`K = ℚ₃(ζ₃)`、`σ_m = τ³`(位数 3、固定体 `ℚ₃(ζ₉)`)。
`i_m = 8`、`e_F = 6`、`p·e_F = 18 = e_L`、`s_m = 16/18 = 8/9`、★`3 ∤ 8`。 -/
theorem not_chain_ledger_of_not_dvd {p i eK : ℕ} (hp : 1 < p) (heK : 0 < eK) (hnd : ¬ p ∣ i)
    {s : ℕ → ℝ} {m : ℕ} (hm : 1 ≤ m)
    (hs : s m = (((p:ℝ) - 1) * (i : ℝ)) / ((p : ℝ) * (eK : ℝ)))
    (hsm : s m ≤ 1)
    (hrec : ∀ j, j < m → (p : ℝ) * s j ≤ s (j + 1)) :
    ¬ (1 - ((p : ℝ) - 2) * s m / ((p : ℝ) - 1) - (∑ j ∈ Finset.range m, s j)
        ≤ ((p : ℝ)⁻¹) ^ m / ((p : ℝ) - 1)) :=
  fun h => hnd (dvd_of_chain_ledger hp heK hm hs hsm hrec h)

def not_chain_ledger_of_not_dvd.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §3 数値の反例(★実在の分岐データ) -/

/-- ★★★実在の分岐データによる反例(★数値版)。

`p = 3`、`L = ℚ₃(ζ₂₇)`、`K = ℚ₃(ζ₃)`、`e_L = 18`、`i(τ) = 2`、`i(τ³) = 8` なので
`s₀ = 2·2/18 = 2/9`、`s₁ = 2·8/18 = 8/9`。★鎖の仮定(`3·s₀ ≤ s₁`、`s₁ ≤ 1`)は全部満たすのに

  総損失の指数 `= 1 − (3−2)·(8/9)/2 − 2/9 = 1/3`  >  `axDecay 3 2` の指数 `3⁻¹/2 = 1/6`。

★★`WildJumpChain.not_ledger_of_ge` は同じデータで `2/9 > 1/6` と測っていた。
★正しい損失項(§測定 1)を入れると `1/3` で、★差はもっと大きい★。 -/
theorem not_chain_ledger_cyclotomic :
    ∃ s : ℕ → ℝ, (∀ j, 0 ≤ s j) ∧ s 1 ≤ 1 ∧ (∀ j, j < 1 → (3:ℝ) * s j ≤ s (j + 1)) ∧
      ¬ (1 - ((3:ℝ) - 2) * s 1 / ((3:ℝ) - 1) - (∑ j ∈ Finset.range 1, s j)
          ≤ ((3:ℝ)⁻¹) ^ 1 / ((3:ℝ) - 1)) := by
  refine ⟨fun j => if j = 0 then 2/9 else 8/9, ?_, ?_, ?_, ?_⟩
  · intro j; dsimp only; split <;> norm_num
  · norm_num
  · intro j hj
    have hj0 : j = 0 := by omega
    subst hj0
    norm_num
  · simp only [Finset.sum_range_one]
    norm_num

def not_chain_ledger_cyclotomic.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §4 ★点 2(`hbase`)—— 木に既に在ったものの橋 -/

section Bridge

variable {L : Type*} [NormedField L]

/-- ★★★点 2(`hbase`)の抽象核。`WildJumpChain.mul_le_of_norm_natCast_le_pow` の★逆★。

`‖(n:L)‖ = ‖π‖^e` と `N·i ≤ e` から `‖(n:L)‖ ≤ (‖π‖^i)^N`。
`n = p`、`N = p−1`、`i = i(σ)`、`e = e_L`、`θ = ‖π‖^i` を入れると

  ★`‖p‖ ≤ θ^{p−1}`  —— ★配られた `hbase` の字面そのもの。

★★入力 `(p−1)i ≤ e_L` は `RamificationJumpBound` が持っている。それどころか
`RamificationJumpBound.norm_natCast_le_pow_of_splits` の結論
`‖(n:L)‖ ≤ ‖π‖^{(n−1)i}` は★既に `hbase` そのもの★で、本補題はその並べ替えである。
⇒ ★★`WildJumpChain` の「木に無い」は★外れ★であった(モジュール docstring の在庫の測定)。 -/
theorem norm_natCast_le_pow_of_mul_le {N i e n : ℕ} {π : L}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1) (he : ‖(n : L)‖ = ‖π‖ ^ e) (h : N * i ≤ e) :
    ‖(n : L)‖ ≤ (‖π‖ ^ i) ^ N := by
  rw [he, ← pow_mul, mul_comm i N]
  exact pow_le_pow_of_le_one (le_of_lt hπ0) hπ1 h

def norm_natCast_le_pow_of_mul_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Bridge

/-! ## §5 ★`θ = 1` での退化(★測定 3 の帰結) -/

/-- ★★`θ = 1` では平均の道の損失は `p` そのもの。

★収縮率の仮定 `∀ z : K.closure, ‖σz − z‖ ≤ θ‖z‖` は `θ < 1` では満たせない
(モジュール docstring §測定 3)ので、この値しか使えない。
⇒ ★`WildDepthFieldDescent.axWildDescent_prime` の定数 `p` に戻る。 -/
theorem orbit_average_loss_at_one {p : ℕ} (hp : 1 ≤ (p:ℝ)) :
    max ((p:ℝ) * (1:ℝ) ^ (p - 2)) 1 = (p:ℝ) := by
  rw [one_pow, mul_one]
  exact max_eq_left hp

def orbit_average_loss_at_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★`θ = 1` では `ε` の縮みが `1`(★稼ぎが無い)。§5 の相方。 -/
theorem contract_gain_at_one {p : ℕ} (hp : 1 ≤ (p:ℝ)) :
    max ((1:ℝ) ^ (p - 1)) ((p:ℝ)⁻¹) = 1 := by
  rw [one_pow]
  refine max_eq_left ?_
  rw [inv_le_one₀ (by linarith)]
  exact hp

def contract_gain_at_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §6 ★点 1 —— 落ちたのは糊だけ(★sharp な形は閉じていない) -/

/-- 抽象核(純 `ℕ`): 最初の跳び `u₁` は指数 `p` の商の跳び `i` 以下なので

  `u₁ ≤ i` かつ `(p−1)i ≤ p·e_K` ⟹ ★`(p−1)·u₁ ≤ p·e_K`。

★★これは★糊だけ★である。中身は `u₁ ≤ i(σ̄)` の方で、そちらは Herbrand
(`φ_{L/K}` が `u ≤ u₁` で恒等)= Serre CL IV §1 Prop. 3 であり、★形式化していない★。 -/
theorem sub_one_mul_first_jump_le {p u i eK : ℕ} (hui : u ≤ i) (h : (p - 1) * i ≤ p * eK) :
    (p - 1) * u ≤ p * eK :=
  le_trans (Nat.mul_le_mul_left _ hui) h

def sub_one_mul_first_jump_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★測定: 「距離だけで Herbrand を代用する」と★どれだけ弱くなるか★(純 `ℕ`)。

ノルムから出るのは正規化された比 `u₁/e_L ≤ i/e_{L^N}`(= `u·e_N ≤ i·e_L`)だけなので、
`(p−1)i ≤ e_N` と合わせても出るのは ★`(p−1)·u₁ ≤ e_L`。
★欲しい `(p−1)u₁ ≤ p·e_K = e_N` より `e_L/e_N = [L:K]/p` 倍だけ弱い。

⇒ ★★「Herbrand は距離では代用できない」ことが機械で測れた形である。 -/
theorem sub_one_mul_first_jump_le_of_metric {p u i eL eN : ℕ} (h0 : 0 < eN)
    (hmet : u * eN ≤ i * eL) (hb : (p - 1) * i ≤ eN) :
    (p - 1) * u ≤ eL := by
  refine Nat.le_of_mul_le_mul_right ?_ h0
  calc (p - 1) * u * eN = (p - 1) * (u * eN) := by ring
    _ ≤ (p - 1) * (i * eL) := Nat.mul_le_mul_left _ hmet
    _ = ((p - 1) * i) * eL := by ring
    _ ≤ eN * eL := Nat.mul_le_mul_right _ hb
    _ = eL * eN := by ring

def sub_one_mul_first_jump_le_of_metric.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §7 ★点 3 —— 足りないのは定数だけである -/

/-- ★★`AxWildDescent K c` は定数について単調(`c ≤ c'` なら `c'` でも成り立つ)。

★★これで点 3 の位置づけがはっきりする: 降下そのものは
`WildDepthFieldDescent.axWildDescent_normInv` が★既に持っている★(深さちょうど 1 減、定数 `p`)。
`AxWildDescent` も `AxTowerDecay.exists_mem_of_descent_budget` も
★狭義の減少しか要求しない★ので「ちょうど 1」も「最初の跳びの層」も要らない。
⇒ ★残っているのは★定数を `axDecay p k` まで絞ること★ただ 1 つで、
★それが §1–§3 のとおり鎖では★できない★。 -/
theorem axWildDescent_mono {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) {c c' : ℕ → ℝ}
    (hle : ∀ k, c k ≤ c' k) (h : AxWildDescent K c) : AxWildDescent K c' := by
  intro ε hε x hxε hdvd
  obtain ⟨x', hlt, hd, hinv⟩ := h ε hε x hxε hdvd
  exact ⟨x', hlt, le_trans hd (mul_le_mul_of_nonneg_right (hle _) hε),
    fun σ => le_trans (hinv σ) (mul_le_mul_of_nonneg_right (hle _) hε)⟩

def axWildDescent_mono.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §8 使っている公理の一覧 -/

#print axioms step_ge
#print axioms sum_le_of_chain_ge
#print axioms chain_ledger_forces_max
#print axioms dvd_of_chain_ledger
#print axioms not_chain_ledger_of_not_dvd
#print axioms not_chain_ledger_cyclotomic
#print axioms norm_natCast_le_pow_of_mul_le
#print axioms orbit_average_loss_at_one
#print axioms contract_gain_at_one
#print axioms sub_one_mul_first_jump_le
#print axioms sub_one_mul_first_jump_le_of_metric
#print axioms axWildDescent_mono

end ABC3.Found.PGC
