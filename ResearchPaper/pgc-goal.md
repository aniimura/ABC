# [pGC] ゴール —— A Version of the Grothendieck Conjecture for p-adic Local Fields

2026-09-05 に作成。`corrhyp-goal.md` / `genell-goal.md` と同じ位置づけ。

## 到達点(原典の項目)

| 節 | 項目 | 状態 |
|---|---|---|
| §1 | Proposition 1.1(円分指標が Γ_K から回復) | sorry |
| §1 | Proposition 1.2(q と [K:ℚ_p] が Γ_K から回復) | sorry・**形式化が偽だった**(第 1012) |
| §1 | Corollary 1.3(惰性群 I_K が Γ_K から決まる) | sorry・**Prop 1.2 へ帰着済み**(第 995) |
| §2 | Proposition 2.1(Γ_K-加群 K̄ の回復) | sorry |
| §2 | Proposition 2.2(O_K̄ と K̄^ の回復) | sorry |
| §2 | Definition 2.3 | **完了** |
| §3 | Corollary 3.1(Hodge-Tate 性の判定) | sorry |
| §3 | Definition 3.2 | **完了** |
| §3 | Corollary 3.3 | sorry |
| §4 | Lemma 4.1 | **完了**(原典で唯一、著者が完全な自己完結証明を与えている) |
| §4 | Theorem 4.2(主定理) | sorry |

★**件数を目標にしない。** 今日いちばん価値があったのは Prop 1.2 の形式化が偽だと
突き止めたこと(7 例目)だが、これは件数を 1 つも進めない。
この木の失敗形は「偽になる」か「自明になる」の二択で、
**退化を見つけることは項目を埋めることと同じ重み**を持つ。

## いま塞いでいるもの(2026-09-05 の測定)

`Skeleton/PGC/Section1.lean`(下流 32 ノード / 12 項目)が pGC 全体の栓である。
§2・§3・§4 はすべてここで止まる(`node tools/frontier.mjs --owner pGC --all`)。

- **Prop 1.1** —— 原典の飛躍は 5 行の補題(階数 1 自由 ℤ_p 加群の同型は指標を決める)。
  ただしその手前の「ℤ_p(1) の同型類が回復できる」に**局所 Tate 双対性**が要り、
  mathlib は**カップ積が無いので述べることすらできない**段階。
  → 退化の検査は済んだ(第 1011、`χ ≢ 1`)。
- **Prop 1.2** —— 下記の経路 C。
- **Cor 1.3** —— Prop 1.2 が付けば自動(第 995 `inertia_recoverable_of_prop12`)。
- **Cor 3.1 / Cor 3.3 / Thm 4.2** —— `K ≠ K'` の witness(Iwasawa 型の同型)が
  無いので反証も届かない。`Check/PGC/RefutationAttempts.lean` の壁。

## ★経路の選択(2026-09-05 に決定)

Prop 1.2 について 3 経路を比較した。

| 経路 | 内容 | 規模 |
|---|---|---|
| A 局所類体論 | `Γ_K^ab ≅ (K^×)^∧` を作る。`K^ab = K_π·K^ur` が本丸 | **既存 Lubin-Tate 塔(70 ファイル)と同規模** |
| B 馴分岐 `I_K/P_K` | Frobenius が q 倍で作用することから読む | ★**Cor 1.3 と循環**(I_K の決定は Prop 1.2 の下流) |
| **C アーベル化** | `N_n(Γ) := #Hom_cont(Γ, ℤ/n)` を不変量に取る | **新規ノード 8〜10 個** |

**C を採る。** 原典は Serre(局所類体論)に投げているが、我々が消費する結論は
「q と [K:ℚ_p] が Γ_K から決まる」だけであり、そこへは初等的に届く。
★逸脱として記録する:原典の論拠(相互律)を経由しない。

### C の主張

**(C-q)** `p ∤ n` のとき `N_n(Γ_K) = n · gcd(n, q−1)`。
ゆえに `q − 1 = max{ n : p ∤ n, N_n(Γ_K) = n² }`。

- 下界は在庫(`exists_surjective_abelianization_zmod_prod_units` の n=1 段)。
- 上界は `K^ur` へ基底変換して Frobenius の余不変量を取る。
  `μ_n ⊂ K^ur`(在庫、第 1005/1006)なので Kummer 理論が使え、
  `𝒪_{K^ur}^×` は n 可除(Hensel)。Frobenius は π に自明・`μ_n` に **q 乗**で作用
  (在庫、第 1010 `exists_frobenius_pow_of_pow_eq_one`)。

**(C-d)** `μ_p ⊆ F` なる `F` で `N_p(Γ_F) = |F^×/(F^×)^p| = p^{[F:ℚ_p]+2}`。
Kummer 双対(全射性は Hilbert 90、mathlib 在庫)+ `unitsSplitEquiv` +
`index_powRange_smallPrincipalUnits`。一般の K へは特性的開部分群
`ker(Γ_K → Γ_K^ab/(p−1))` に落として戻す(部分群対応は在庫、第 1000 前後)。
**Cor 1.3 を使わないので循環しない。**

### ★不変量の定義(2026-09-05 決定)

`N_n(Γ) := Nat.card {f : Γ →* Multiplicative (ZMod n) // IsOpen (ker f)}`。

- ✗ `ContinuousMonoidHom` —— **mathlib に `TopologicalSpace (ZMod n)` が無い**(実測)
- ✗ `Abelianization Γ` の `(powMonoidHom n).range.index` —— **偽になりうる**。
  `Abelianization` は代数的交換子群による商で位相的閉包を取らない。
  副有限群では有限指数部分群が開とは限らない
- ✓ `IsOpen (ker f)` —— 離散位相を持ち込まずに済み、移送は α が開写像であることだけ

### C のノード(★`Found/PGC/` に置く。原典の項目ではないので `.src` を持てない——D1)

| | ファイル | 内容 | 見積 | 依存 |
|---|---|---|---|---|
| **A** | `ContinuousHomCount.lean` | `contHom` / `contHomCard` / `contHomCard_congr` | 小 | — |
| **B** | `UnitsPowPrimeToP.lean` | `p∤n` で `[𝒪^×:(𝒪^×)^n] = gcd(n,q−1)`、`[K^×:(K^×)^n] = n·gcd` | 中 | — |
| **C** | `KummerDuality.lean` | `contHomCard Γ_F n = [F^×:(F^×)^n]`(★**要石**) | 大 | — |
| **D** | `HerbrandIndex.lean` + `UnitsPowP.lean` | 有限指数での Herbrand 商の不変性、`[𝒪^×:(𝒪^×)^p] = p^d·#μ_p` | 中〜大 | — |
| **E** | `UnramifiedClosureRoots.lean` | `K^ur` のノルム 1 の元は `p∤n` 乗根を持つ | 中〜大 | — |
| **F** | `InertiaKummer.lean` | 惰性群からの連続準同型、Frobenius の余不変量 | 大 | C, E |
| **G** | `ResidueCardTransport.lean` | (C-q) 総組み立て、`residueCard_eq_of_absGal_equiv` | 中 | A,B,C,E,F |
| **H** | `DegreeTransport.lean` | (C-d) 本体、`finrank_eq_of_absGal_equiv` | 中 | A,C,D |
| **I** | `Interface` の修理 | `ResidueCardinality.card_congr`(第 1012、7 例目) | 中 | — |
| **J** | Skeleton 接続 | Prop 1.2 を `⟨G,H⟩` で埋める。★**Cor 1.3 も同時に落ちる** | 小 | G,H,I |

★**A〜E と I は互いに依存しない**ので同時に配れる。

### ★落とし穴 2 つと、その回避(2026-09-05 実測)

1. **Hilbert 90 は `[FiniteDimensional K L]` 必須**で `Γ_F` に当たらない
   (mathlib 自身が無限次への拡張を TODO にしている)。
   → **使わない。** `exists_root_adjoin_eq_top_of_isCyclic` を、各 `f` が切る
   **有限次**巡回拡大に当てる。Γ_F は無限次のままでよい。
2. **`X^n − C a` の既約判定が奇数 `n` 限定**(`TODO: criteria for even n`)。
   → **既約性をどこでも使わない。** Kummer 指標を素手で定義する
   (`κ_a(σ) := σ(α)/α`、α は任意の根)。**これで `p = 2` も素通しできる。**
   ★`autEquivZmod` / `autEquivRootsOfUnity` を使う設計にすると偶数の壁に当たる。

### ★残るリスク

**F1**(`Γ_{K^ur} ≃ₜ* (unramifiedClosure K).fixingSubgroup`)。在庫の
`fixingSubgroupContinuousMulEquiv`(第 995)は有限次限定で、`K^ur` は無限次。
両方向とも作り直しになる(順方向は塔、逆方向は `finiteDimensional_sup` が
有限性を要求する)。**波 2 の前に小ノードで測ること**(`decisions-pending.md` D3)。

## 済んだ土台(2026-09-05 までに sorry 無しで構成)

不分岐拡大理論(存在・一意性・`Gal ≅ ℤ/n`・`K^ur`・`Γ_K/I_K ≅ Gal(K^ur/K)`)/
Lubin-Tate 主定理 `Gal(K(Λ_n)/K) ≅ (𝒪_K/π^n)^×` と次数・完全分岐 /
`K(Λ_n) ⊓ K^ur = K` / 有限段の同時全射 `Γ_K^{ab} ↠ ℤ/m × (𝒪_K/π^n)^×` /
`μ_m ⊆ K^ur` / Frobenius は `μ_m` に q 乗で作用 / `Γ_{L_H} ≃ₜ* H` /
Γ_{ℚ_p} は非可換 / p 進対数(準同型性・単射性・全射性)

## 進め方

`ResearchPaper/orchestration.md` に従う。`node tools/frontier.mjs --owner pGC`
が次の持ち場を出し、`node tools/brief.mjs --node <rel>` が agent に渡す形にする。
同時起動は 5 まで。全体ゲートは実装 agent が全員止まってから。
