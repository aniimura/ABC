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

**How to apply:**
- ★**件数は `node tools/check.mjs` の「Interface 実装待ち」が数える。**
  2026-08-17 時点で 26 件(Arakelov 9 + Galois 8 + 既存 9)。
  ★★増えて見えるのは**後退ではない**——畳まれていたものが数えられるようになっただけ。
- ★★**posit は最小にする。** `E[n]` は mathlib の `W.toAffine.Point`
  (`AddCommGroup`、★`[DecidableEq K]` が要る)から**今すぐ書ける**ので
  `torsionPoints` は `def` にし、**構造定理だけ** posit した。同じ判断を他でもすること。
- ★★★**`ℤ_[l]` は `[Fact l.Prime]` を要求する**——構造体のフィールドでも束縛子が要る。
- 詳細は `ResearchPaper/genell-goal.md` §9-20。
- 関連: [[genell-track-b]] / [[lean-build-check-discipline]] /
  [[parallel-session-sweeps-my-files]]
