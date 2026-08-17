---
name: arakelov-galois-skeleton-counts
description: Arakelov 理論は 9 件、Galois 表現は 8 件に割ってある。律速は B1 層のテンソル積・C2 ℙⁿ の点の関手・G1 E[n]≅(ℤ/n)²。件数は check.mjs の「Interface 実装待ち」が数える。
metadata:
  type: project
---

**S2・S3 の基礎理論は `Interface/Arakelov/` と `Interface/GaloisRep/` に割ってある**
(2026-08-17 作成)。★それ以前は `Interface/GenEll/` の `waiting` **大きな文字列 1 本**で、
「何本作れば埋まるか」が数えられなかった。

## ★★★Arakelov 理論 = 9 件(`Interface/Arakelov/`)

| # | obligation | 状態 |
|---|---|---|
| B1 | `Pic(X)` 可逆層の群 | ★★★**律速**。`IsLocallyFree` は有るが**階数**と**層のテンソル積**が無い(前層版は有る→層化で届く) |
| B2 | `𝒪_X(D)` Cartier → 可逆層 | (B1) に従属。`comap_mul` は我々が証明済 |
| B3 | `Pic(Spec 𝓞_F) ≅ ClassGroup` | (B1) が入れば機械的 |
| C1 | `X^arc` の位相・`ι_X`・コンパクト | ★★**構成済**(`Found/GenEll/ArcModel.lean`) |
| C2 | ℤ-固有 ⇒ 射影埋め込み | ★★★**律速**。`ℙⁿ` の**点の関手**が無い |
| C3 | 解析化と hermitian 計量 | (B1) に従属。`IsConjInvariant` は定式化済 |
| D1 | `APic(X)` | (B1)+(C3) に従属 |
| D2 | `APic(Spec 𝓞_F)` と `deg_F` | `ADiv`/`deg_F`/`APrc` 実装済、底変換不変性も証明済。橋だけ |
| D3 | `ht_L̄` と `Prop 1.4` | ★★★**`U_X(ℚ̄)` では構成済**。(B2) が入れば全域化 |

## ★★★Galois 表現 = 8 件(`Interface/GaloisRep/`)

| # | obligation | 状態 |
|---|---|---|
| G1 | `E[n] ≅ (ℤ/n)²` | ★★★**S3 の入口かつ壁**。mathlib にも **FLT にも無い**(FLT は sorry) |
| G2 | Tate 加群 `T_l E ≅ ℤ_l²` | (G1) が入れば機械的 |
| G3 | `ρ_{E,l} : Gal → GL₂(ℤ_l)` | 行き先と定義域は書ける。★**Weil 対**が無い |
| G4 | `mod l` 表現 | `PadicInt.toZMod` は有る。★`Lemma 3.1` は実装済 |
| G5 | 像が `SL₂` を含む(`Theorem 3.8`) | S3 の最後の段 |
| G6 | Tate 曲線と局所高さ(`Definition 3.3`) | 両方に無い。典拠 [FC] III, Cor 7.3 |
| G7 | 半安定還元と `𝓞_L` 上のモデル | `Reduction.lean` は有るが Néron モデルは無い |
| G8 | Faltings 高さ | ★★★**Arakelov 側との合流点**((D2) が要る) |

**Why:** ★★★**`E[n] ≅ (ℤ/n)²` が入口である**——これが無いと `GL₂` の **`2`** が書けない。
★Arakelov 側は **B1(層のテンソル積)と C2(`ℙⁿ` の点の関手)の 2 本**が律速で、
それ以外はすべて従属か既に構成済みである。


## ★★★2026-08-17 追記: 17 件すべてを**退化不能**にした

★★★**最初のスケルトンは「中身ゼロの witness」で埋まった**(実測):

| obligation | 退化 witness |
|---|---|
| B1 | `Pic := PUnit`(自明群) |
| B2 | `𝒪(D) := 1` |
| C1 | `logMetric`/位相を離散に |
| C3 | `logMetric := 0` |
| D1 | `APic := PUnit` |
| D2・D3 | `degF := 0` / `height := 0` |
| G3 | `rep := 1`(自明表現) |
| G6・G7・G8 | `LocalField := PUnit` / `SemiStable := True` / `htFalt = degInf = 0` |

★★**このままでは「9/9 達成」が数学をせずに書けてしまう。**全件を塞いだ。

### ★塞ぐ 2 つの手

1. ★★★**posit を mathlib の対象に接地する**——自前の型は 1 点に潰せるが、
   mathlib の型は潰せない。G6 を `HasSplitMultiplicativeReduction` に、
   B1 を `CommRing.Pic` に置き換えた(**どちらも mathlib に在った**)。
2. ★**値を別の場に縛る公理を足す**——
   `localHeight = v(q)` / `omega` の階数 = 1 / `logMetric(scale c m) = logMetric m + c` /
   `forgetMetric (ofMetric L m) = L` / `degF(scale c m) = degF m − c` /
   `deg∞ ≥ 局所高さ·log2` / `det ∘ rep` が全射。

### ★★★必ず負の対照を取る

退化 witness を**実際に書いて落ちることを確認**した。
★例: `localHeight := 1` → `⊢ False`、`omega := PUnit` → 階数 0 で矛盾、
`logMetric := 0` → `⊢ 0 = c`、`Pic := PUnit` → `Unique (CommRing.Pic R)` が無い。

## ★★★★★C1 は達成した(2026-08-17)——Arakelov **1/9**

`ArcSpaceData.nonvacuous` が取れた(`Found/Arakelov/` 14 ファイル、sorry 0、標準 3 公理のみ)。
★`check.mjs` の「Interface 実装待ち」は **26 → 25 件**。

### ★3 段で落ちた

| 段 | 定理 |
|---|---|
| A | `arcTopology_opens_of_affine`——`Spec A` の**任意の**開部分スキーム |
| B 前半 | `isOpenMap_comp_of_isAffine`——任意のアフィン標的 |
| B 後半 | `arcTopology_openImmersion`——一般の開埋め込み |

★★★段 B の核心は「**`U ⊓ O` を経由する**」ことだった:

    (· ≫ U.ι) ⁻¹' ((· ≫ O.ι) '' V) = (· ≫ homOfLE) '' ((· ≫ homOfLE) ⁻¹' V)

`arcTopology = ⨆`(アフィン chart)なので開性は chart ごとに降ろせ(`isOpen_iSup_iff`)、
各 chart は**アフィン**だから段 A が効く。

### ★★★GAGA は要らなかった

当初「複素解析空間と GAGA が要る」と見積もったのは**誤り**だった。
実際に要ったのは**商位相と多項式の連続性**だけである。

## ★★★★退化封じの見落とし——「型を固定する」だけでは足りない

★★★C1 は当初、次の退化 witness を**通してしまっていた**:

    evalAffine := fun _ _ _ => 0
    topology   := fun _ => ⊤        -- 密着位相

理由: `induced (fun _ _ => 0) Pi = ⊤` なので `topology_affine` は
「`topology = ⊤`」を要求するだけになり、`topology_openImmersion` も
`induced g ⊤ = ⊤` で自明に成り立つ。

★`equivComplexPoints` は**台の型**を固定するが、**`evalAffine` の値**は固定しない。
★塞ぎ方: `evalAffine_spec`(評価は `Spec.preimage` が与える環準同型)を足した。
実装側では **`rfl`** で通る。

★★★**教訓: 退化封じでは「その構造で値が自由に選べるフィールド」を列挙すること。**
型を固定しただけで満足しない。

## ★★★Lean 実装で 4 度かかった罠

1. ★★**`simp only [Function.comp_def]` は `Scheme.Hom.mk` まで展開して壊す**。
   `rw [h]` を使い、**`h` を「合成の形」で述べる**(ゴールが既に合成形)。
2. 位相のインスタンスは `letI` で**明示的に**入れる(`⨆` の成分は自動で決まらない)。
3. ★★★**`TopologicalSpace` の `≤` は「細かい」**。`le_def` の表示は
   `IsOpen ≤ IsOpen` で**左右が読めない**——**`⊥` が離散かを試して**確定させる。
   ★一度これで `⨅`/`⨆` を逆にした。
4. `Spec ℂ` の点の型は **`↥(Spec (CommRingCat.of ℂ))`** で
   `PrimeSpectrum (…)` と**構文的に別**(defeq だが `rw` は噛まない)。

## ★★★★★B1 は 14 ブロックまで積んだ(2026-08-17)——`CommGroup` の 5 公理が揃った

`Found/Arakelov/` に 7 ファイル(`PicPresheafTensor` … `PicType`)、すべて sorry 0。

| 公理 | 定理 |
|---|---|
| 乗法 | `tensorModules`(前層でテンソル → 層化) |
| 可換 | `tensorModulesComm` |
| 単位元 | `tensorUnitLeft` / `InvertibleSheaf.one` |
| 結合律 | `tensorModulesAssoc` |
| 逆元 | `InvertibleSheaf.symm` |

★★★**残り 2 点**: (a) テンソル積が可逆層で閉じること
(`IsLocallyRankOne` を「`Over V` 上の前層同型」へ強める必要がある)、
(b) 同型類への商と `CommGroup` インスタンス(`Shrink` の本質的小ささも)。

### ★★★見積りを 4 度外した(すべて「無い」と決めてから在庫が見つかった)

| 何 | 当初の見積り | 実際 |
|---|---|---|
| C1 の `X^arc` | 複素解析空間と GAGA | ★商位相と多項式の連続性だけ |
| 開集合への制限 | 自作 | ★`SheafOfModules.over` が在った |
| 制限とテンソルの両立 | 自作 | ★`pushforward₀OfCommRingCat.Monoidal` が在った |
| 結合律 | 内部 Hom が要る(mathlib PR 級) | ★局所論法で回避、`Over` すら不要だった |

★★★**「無い」と決める前に測る。**

### ★★★★★Lean の罠(この turn で 5 度)

★★**インスタンス束縛子は「型の書き方の違い」をまたげない。**
`X.PresheafOfModules` と `PresheafOfModules (X.presheaf ⋙ forget₂ _ _)` は
`rfl` で等しいが、インスタンス探索は片方でしか成功しない。
★**対処: 条件を `[..]` でなく `(h : ..)` で受ける。**
★別名(`abbrev` でも)を経由すると届かないことがある——**素の形で書く**。
★依存を書き換えたら `lake build` で先に olean を作り直す
(さもないと名前付き引数が古い署名で解決される)。

## ★★★★★★B1 は 17 ブロック——`CommGroup (PicPre X)` まで到達(2026-08-17)

`Found/Arakelov/` は **28 ファイル・251 宣言・sorry 0**。

### ★★★★★決定打は「**層化を通さない**」だった

第 13 ブロックまで「前層でテンソル → 層化」で組んでいたが、
★層化を通すたびに「局所自明性が保たれるか」を示す必要があった。

★★★**第 16 ブロックで前層の段に移したところ、群法則が無料になった**:

| 何 | 前層の段 | 層化を通す段 |
|---|---|---|
| モノイダル構造 | ★mathlib が持っている | 無い(13 ブロックかけて作った) |
| 結合律・単位律・可換律 | ★★**`α_` / `λ_` / `ρ_` / `β_` そのもの** | 局所論法が要る |

★局所自明な前層は自動的に層(層条件は局所的)なので数学的に一致する。

### ★★B1 の残り 3 点(すべて実測済み)

1. 引き戻し `f^*` —— mathlib は `Scheme.Modules.pullback` を持つが
   ★★**モノイダル性は未登録**。我々の `Pic` は前層側なので橋も要る
2. **局所自明な前層は層**(`Interface` の `sheafOf` に要る)
   ——「層条件は局所的」を site の言葉で。★mathlib に該当補題なし
3. `equivPicRing` —— `Tilde.lean` の `fromTildeΓ` が材料

★★★2 は「層化が局所自明性を保つ」と同値で、**`sheafOf` を層化で定義しても回避できない**
(`sheafifyTensorLeft/Right` の仮定が層化後の対象に掛かるため)。

### ★★★★★見積りを 5 度外した

| 何 | 当初 | 実際 |
|---|---|---|
| C1 の `X^arc` | 複素解析空間と GAGA | ★商位相と多項式の連続性だけ |
| 開集合への制限 | 自作 | ★`SheafOfModules.over` が在った |
| 制限とテンソルの両立 | 自作 | ★`pushforward₀OfCommRingCat.Monoidal` が在った |
| 結合律 | 内部 Hom(mathlib PR 級) | ★局所論法で回避、`Over` すら不要 |
| 群構造 | 層化を通す必要がある | ★★**前層の段で組めば無料** |

★★★**「無い」と決める前に測る。**

## ★★★★★★★2026-08-17 の訂正: 「局所自明な前層は層」は**偽**だった

★前層の段で `Pic` を組もうとしたが、根拠が偽だった。反例:

    (𝒪(1) ⊗_pre 𝒪(-1))(ℙ¹) = H⁰(𝒪(1)) ⊗_k H⁰(𝒪(-1)) = k² ⊗ 0 = 0  ≠  k

★★前層テンソルは**各開集合ごと**なので、`F ⊗_pre G ≅ 𝟙_` は
**本物の直線束では成り立たない**。★「層条件は局所的」も一般には偽
(貼り合わせを示すのに貼り合わせを使う循環)。

★★★**4 ファイル(`PicGroup` / `PicQuotient` / `PicInterface` / `PicAffine`)を撤回した。**
第 1–15 ブロック(層側)は無傷。

### ★★★★★教訓——両方向に測る

| 姿勢 | 実績 |
|---|---|
| 「**無い**」と決める前に測る | ★7 度成功(GAGA / 制限 / 制限とテンソル / 内部 Hom / 層化 / 局所全単射の保存) |
| 「**有る(自明)**」と決める前に測る | ★★★**1 度失敗**——具体例で検算していれば 5 ブロック前に気づけた |

★★★★**「自動的に成り立つ」と思った補題こそ、具体例で検算する。**

## ★★B1 の最後の壁と、その部品(2026-08-17 実測)

壁: **層化は局所自明性を保つ**。★部品は**すべて mathlib に在った**:

| 部品 | 在庫 |
|---|---|
| `Over.forget` が cocontinuous | `Sites/Over.lean` |
| 制限が局所全射/単射を保つ | ★`Sites/PreservesLocallyBijective.lean` |
| 層どうしの局所全単射は同型 | ★`Sites/LocallyBijective.lean` の `Sheaf.isLocallyBijective_iff_isIso` |

★残るのは `PresheafOfModules` の射を `Sheaf J A` の射として包み直す手間だけ。

## ★★★★★★★★B1 の中核が完成した(2026-08-17)——`CommGroup (PicSheaf X)`

`Found/Arakelov/` 15 ファイル、すべて sorry 0。★★**層の側**で組んである
(前層の段で組もうとした最初の試みは誤りで撤回済み)。

| 段 | 定理 |
|---|---|
| 壁 | ★`isLocallyTrivial_sheafify`(層化は局所自明性を保つ) |
| 閉性 | ★`isLocallyTrivial_tensorModules` |
| 並べ替え | ★`tensorRearrange` |
| 群 | ★★★`instance : CommGroup (PicSheaf X)` |

### ★★壁の破り方——4 部品すべて mathlib に在った

    η_P|_V は局所全単射(`Sites/PreservesLocallyBijective.lean`)
    両側とも層(`Functor.op_comp_isSheaf`)
    ⟹ 同型(`Sheaf.isLocallyBijective_iff_isIso`)
    ⟹ 前層加群の射としても同型(`toPresheaf` が同型を反映)

★★★**「mathlib に該当補題なし」と 2 度判定したが、いずれも誤りだった。**

### ★B1 の残り 2 点——どちらも「関手がテンソル積を保つ」

- 引き戻し `f^*`: `Scheme.Modules.pullback` は在るが**モノイダル性が未登録**
- `equivPicRing`: `tilde` の**モノイダル性**

★mathlib は `pushforward₀OfCommRingCat` について既に作っているので同じ型の仕事。

## ★★グラフにも 28 件が出るようになった(2026-08-17)

`dependency-graph.html` の節点語彙は**原文の項目番号**なので、
`GenEll Definition 1.1` に 8 本ぶら下がっていても **1 節点に潰れて見えなかった**。
★`tools/graph-html.mjs` に `Interface` タグの義務節点を足した(28 件、埋まった 3)。
★★被覆率(節点 1015 / 着地 10 / 張った 54)は**汚していない**——義務は統計から除く。

**How to apply:**
- ★**件数は `node tools/check.mjs` の「Interface 実装待ち」が数える。**
  2026-08-17 時点で **25 件**(C1 達成で 26 から 1 減った)。
  ★`node tools/graph-html.mjs` の「Interface の義務: 3 / 28」も同じ単位である
  (28 = 25 待ち + 3 埋まった)。
  ★★増えて見えるのは**後退ではない**——畳まれていたものが数えられるようになっただけ。
- ★★**posit は最小にする。** `E[n]` は mathlib の `W.toAffine.Point`
  (`AddCommGroup`、★`[DecidableEq K]` が要る)から**今すぐ書ける**ので
  `torsionPoints` は `def` にし、**構造定理だけ** posit した。同じ判断を他でもすること。
- ★★★**`ℤ_[l]` は `[Fact l.Prime]` を要求する**——構造体のフィールドでも束縛子が要る。
- 詳細は `ResearchPaper/genell-goal.md` §9-20 / §9-23 / §9-24。
- ★★★**退化封じは実装の前に必ずやる**——後だと「埋めた」ものが無内容になる。
- ★★**負の対照(退化 witness を書いて落ちることの確認)を省かない。**
- 関連: [[genell-track-b]] / [[lean-build-check-discipline]] /
  [[parallel-session-sweeps-my-files]]
