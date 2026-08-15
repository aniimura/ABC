# 層0の仕分け——いま着手できるもの(2026-08-15)

依存グラフの層0(= 記録された辺を持たない最左層)を、壁の種類で仕分けた。
目的は**工数の把握**であり網羅ではない。

## ★壁の種類は5つになった

| | 種類 | 実例 | 重さ |
|---|---|---|---|
| 1 | **mathlib の不在** | p進対数(評価 API がノルム体に原理的に不適用) | 作業量 |
| 2 | **Lean の透明度** | tempered 群、packet-normalization | 作業量 |
| 3 | **どれでもない(安い)** | 対数体積、型の訂正。見分け方: 原文が「elementary」と書く / 原文に当たれば決まる | 小 |
| 4 | **語彙の不在** | procession-normalization | ★最も重い。中層を作るまで手が付かない |
| 5 | ★**mathlib にそのものがある** | `integral` / `saturated` / `essential image` | ★最小 |

## 層0 = 129 節点(§0語 47 / 概念 18 / 番号付き項目 64)

### §0 語 47(実語 38 + 語未同定のプレースホルダ 9)

| 種類 | 件数 | 備考 |
|---:|---:|---|
| 5. mathlib にそのものがある | **3** | 下表。原文の定義文と突き合わせ済み |
| 3. 安い | 1 | `sharp`(1行で書ける) |
| 4. 語彙の不在(grep 0件) | 22 | scheme-like / slim / sturdy / nexus / id-rigid / isomorph / fsm-morphism … |
| 未分類 | 12 | 当たりは出たが概念が別物の疑い(`order` 2347件、`irreducible` 644件 …) |

**種類5の3件(原文の定義と突き合わせ済み)**:

| 原文 | mathlib |
|---|---|
| [FrdI] §0 `integral`: `M → M^gp` が単射 | `GrothendieckGroup.of` + `of_injective [IsCancelMul M]` |
| [FrdI] §0 `saturated`: `n·a` が像に入れば `a` も像に入る | `Submonoid.PowSaturated`(`GroupTheory/Subgroup/Saturated.lean:31`) |
| [EtTh] §0 `essential image` | `CategoryTheory.Functor.essImage` |

★**`sharp` は「語の一致は当てにならない」の実例。** grep は 17 件出したが**全部「sharp bound」**(Stirling の評価など)で、[FrdI] §0 の `M^± = 0` とは別物。mathlib に名前は無いが `∀ a, IsUnit a → a = 1` の**1行**なので種類3。

★**mathlib は既に `IsAffineMonoid = IsCancelMul + Monoid.FG + IsMulTorsionFree` として束ねている**
(`Algebra/AffineMonoid/Basic.lean`)。**[FrdI] §0 の語彙の土台が既にある。**

### 概念 18

| 種類 | 件数 | 内訳 |
|---:|---:|---|
| 着手済み | **1** | `log-shell`(ℚ_p について `Found/IUTchIII/PadicLog.lean`、`sorry` 無し) |
| 3. 安い | 2 | `holomorphic hull`(p.127 で `λ·𝒪` と確定、目視済み)/ `field of moduli`(p.113、目視済み) |
| 未分類 / 1 | 2 | `stable reduction`(隣接物 `IsGoodReduction` あり)/ `Galois category`(`FintypeCat` が型に焼き込まれている) |
| 4. 語彙の不在 | 13 | cyclotome / hyperbolic orbicurve / semi-graph / anabelioid / semi-graph of anabelioids / Frobenioid / model Frobenioid / pre-Frobenioid / prime-strip / Hodge theater / Θ-link / mono-theta environment / multiradial |

### 番号付き項目 64

**★親の独立測定(全数、機械)**:

| | 件数 | 意味 |
|---|---:|---|
| `[cf.` を1つも含まない | **41** | **真の葉** |
| `[cf.` を含む | 16 | 引いている先は**外部の古典文献**(`[FC]` Faltings–Chai / `[SGA2]` / `[Serre2]` / `[Serre3]` / `[Ill]` Illusie)、または番号のない論文参照(`[Mzk14]` 等) |
| 解決できず | 6 | 宣言の照合に失敗 |

**子の標本測定(6件、9.4%)**: 着手できるもの 1件(GenEll Lemma 3.1)。→ **64件中 10件前後と推定**(★n=6 なので幅は大きい)。

## ★★いま着手できるもの —— 6件

1. **[FrdI] §0 `integral`(モノイド)** —— `IsCancelMul` + `GrothendieckGroup.of_injective`。原文と一致確認済み
2. **[FrdI] §0 `saturated`(モノイド)** —— `Submonoid.PowSaturated`。原文と一致確認済み
3. **[FrdI] §0 `sharp`(モノイド)** —— `∀ a, IsUnit a → a = 1` の1行。★[FrdI] §0 の中心語
4. **[EtTh] §0 `essential image`** —— `Functor.essImage` そのもの
5. **[GenEll] Lemma 3.1(SL₂ の構造)** —— 純古典。IUT 語彙ゼロ
6. **[IUTchIII] Remark 3.9.5 (i) の `holomorphic hull`** —— p.127 で `λ·𝒪` と確定済み(目視)

## ★★層0は「完成」ではない。理由が重要である

**「依存が無い」は達成されている**(層0の定義そのもの)。**「着手できる」は部分的**である。

> ★**層0であることは、易しいことを意味しない。**
> [EtTh] Proposition 2.2 は層0にありながら `Δ_X`(幾何的基本群)を要る。
> `Δ_X` は §0 でも番号付き項目でもないので辺にならない——**それだけの理由で層0にいる。**

**層0は「どこから始められるか」の下界であって、難易度の順位ではない。**
§0 をグラフに入れて `[FrdI] Definition 1.1` が層0→層1に移ったのと同じ現象が、まだ残っている。
次に埋めるべきは `Δ_X`・`GFG-type` のような**論文固有の記号定義**である。

## ★次に着手すべき1件

**[FrdI] §0 の3語(`integral` / `saturated` / `sharp`)。**

1. 3語とも原文の定義を確認済みで、mathlib に2つは**そのもの**、1つは**1行**。種類5と種類3しかない
2. mathlib が既に `IsAffineMonoid` として束ねている。**土台がある**
3. ★★**これが Frobenioid の入口である。** `Frobenioid` は `[IUTchIII] Definition 3.8` の解消条件であり、
   [FrdI] は機械概算で**最も密**(到達45項目・辺184本・平均出次数 4.1)

> ★**層0の最も安い部分が、中層の最大の壁の入口になっている。**
> これまで「中層は手が付かない」と報告してきたが、**その最下段は mathlib にある。**

## ★残る穴の規模を測った——本文語彙は §0 の 0.2 倍(2026-08-15、親)

子が「次に埋めるべきは `Δ_X`・`GFG-type` のような**論文固有の記号定義**」と指摘した。
その規模を測った(§0 でも番号付き項目でもない場所で定義文型に当たる語):

| | 語数 |
|---|---:|
| §0 語 | **174** |
| **本文語**(§0 でも番号付き項目でもない場所) | **39** |

論文別の本文語(上位): NodNon 13 / CombGC 7 / IUTchIII 4 / IUTchII 3(coric・uniradial・**multiradial**)/
AbsTopIII 3 / Config 3 / MT 3 / IUTchI 2(**conjugate synchronization**・log-theta-lattices)/ FrdI 1(pre-divisorial)

**★判断: 本文語彙の層は今は入れない。**
- 規模が §0 の **0.2 倍**しかない。§0 を入れたときのような効果は見込めない
- ★しかも**この 39 は下界**である。`GFG-type` は
  「… as a profinite group of [almost pro-Σ] **GFG-type** [where …]」という
  **角括弧の挿入を含む**形で定義されており、我々の抽出パターンが捕まえていない
  (AbsTopI の本文語が 0 と出ているのはそのため)
- 入れるなら抽出の精度を上げる必要があり、**その工数に見合う効果が測れていない**

★G7 の基準を当てた: 効果が測れていないものに制度(この場合は抽出層)を足さない。
**代わりにこの測定を残す。** 必要になったら `tools/section0.mjs` の
`section0Terms` を全文に広げれば足りる(§0 の行範囲判定を外すだけ)。
