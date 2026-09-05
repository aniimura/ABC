---
name: pdftotext-two-implementations-hazard
description: pdftotext が 2 実装ある。★罠はシェルではなくキャッシュ。第 1024 で check.mjs 側が対策済み
metadata:
  type: reference
---

この機械には `pdftotext` が 2 つある(2026-09-05 実測):

| | 実装 | 場所 |
|---|---|---|
| Git Bash の PATH 先頭 | **Xpdf 4.00**(較正済み) | `C:\Program Files\Git\mingw64\bin\pdftotext.exe` |
| PowerShell の PATH 先頭 | poppler 25.07.0 | WinGet の Poppler |

**2,157 頁のうち 1,718 頁(79.6%)で違うテキストが出る。**
最多は `´etale`(Xpdf・分解)対 `étale`(poppler・合成)。
`Ŝ` は poppler では `U+0002`+`S`(制御文字が本文に入る)。
NG が 163 件しか鳴らなかったのは**逐語照合が原文のごく一部しか見ていない**からで、
**照合を増やすほど地雷は大きくなる**。

## ★本当の罠はシェルではなくキャッシュだった(2026-09-05 の訂正)

当初この欄には「Git Bash から走らせること」と書いていた。**不十分だった。**
`.cache/pdf-pages.json` の鍵が `パス#頁#mtime#size` で**実装を持っていなかった**ため:

| `pdftotext` | 頁キャッシュ | NG | 時間 |
|---|---|---|---|
| Xpdf(**正しい**) | poppler 産 | **175** | **2 秒** |

**正しいシェルで正しい実装を使っていても 175 件**になり、しかも 2 秒で返るので
「速い＝正常」に見える。**シェルを直しても直らない。**

## 第 1024 で `check.mjs` 側が対策済み

1. `pdftotext -v` から実装を**同定**する
   (★poppler は `-v` を **stderr** に出して**終了コード 0** で返るので `execFileSync` では
   同定できない。`spawnSync` が要る。selftest D44 で固定)
2. `PDFTOTEXT_CALIBRATED = 'Xpdf 4.00'` を、PATH 全体 **+ PATH に無い既知の設置場所**
   (`%ProgramFiles%\Git\mingw64\bin`)から**版で選ぶ**(順序では選ばない)
3. **キャッシュの鍵に実装を入れる**(`{ self, pdftotext, pages }`)

→ PowerShell から走らせても NG 14 のまま。poppler 産キャッシュが残っていても作り直す。

**How to apply:** `node tools/check.mjs --pdftotext` でどれを使うか見られる(0.17 秒)。
較正済みが無ければ stderr に警告が出る(止まりはしない)。
`ABC3_PDFTOTEXT` で明示指定でき、指定が壊れていても PATH へは落ちない。

## ★同じ穴が `0_Source/*.txt` にも開いている(未対策、backlog M11)

`check.mjs` は PDF を直読するが、**他のツールは `.txt` を読む**(137 本)。
`hedge-index` / `paper-items` / `bibmap` / `cycle-probe` / `full-graph` /
`frdi-progress` / `genell-progress` が消費し、★`hedge-index` は
**CLAUDE.md が「着手前に必ず数える」と定めている道具**。
アクセントの指紋で仕分けると **Xpdf 風 111 / poppler 風 3 / 判定不能 23** ——
**既に混ざっている**。`.txt` は人が手で作るので第 1024 の修理の外。

## ★おまけ: `find` は junction を降りない

`.txt` を 0 件と答えるが実際は 137 本ある。**「無い」を `find` で判定しない**
——[[mathlib-cohomology-inventory-2026-09-05]] と同じ「不在の誤り」を道具側で踏む。
